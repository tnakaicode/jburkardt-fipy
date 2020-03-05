#! /usr/bin/env python
#
from fenics import *

def step7 ( my_grid, my_degree, my_case ):

#*****************************************************************************80
#
## step7 compares error norm and DPG error indicators for a Poisson problem.
#
#  Discussion:
#
#    Use the mixed DPG method to solve a Poisson equation:
#
#      - div kappa grad u = f       in Omega = [0,1]x[0,1]
#                       u = u_exact on dOmega
#
#    where:
#
#      kappa(x,y) = 1
#
#    if my_case == 1:
#
#      u_exact(x,y) = x*(1-x)*y*(1-y)
#      f(x,y)       = 2*y*(1-y)+2*x*(1-x)
#
#    if my_case == 2:
#
#      u_exact(x,y) = exp(x^2+y^2)*sin(x+y)
#      f(x,y)       = -2*exp(x^2+y^2)*(sin(x+y)*(1+2*(x^2+y^2)+2*cos(x+y)*(x+y))
#
#    if my_case == 3:
#
#      u_exact(x,y) = 2 * (1+y) / ( (3+x)^2 + (1+y)^2 )
#      f(x,y)       = 0
#
#    To my surprise, it seems that the degree of the Q space and element
#    must be set to 1.  Using a value of 2 results in nan values.
#    JVB, 13 November 2018.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    13 November 2018
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer MY_GRID, the number of (pairs of) elements to use in the 
#    X and Y directions.  Number of triangular elements is 2*my_grid*my_grid.
#
#    Input, integer MY_DEGREE, specifies approximation degrees.
#    Error indicators: MY_DEGREE+2
#    Primary variable: MY_DEGREE+1
#    Fluxes:           1
#
#    Input, integer MY_CASE, chooses the exact solution.
#    1 <= MY_CASE <= 3.
#
  import matplotlib.pyplot as plt
#
#  Report input.
#
  print ( '' )
  print ( '  Exact case #%d' % ( my_case ) )
  print ( '  %dx%d mesh on the unit square' % ( my_grid, my_grid ) )
  print ( '  Approximation degree %d:' % ( my_degree ) )
#
#  Mesh the unit square.
#
  my_mesh = UnitSquareMesh ( my_grid, my_grid, 'right/left' )
#
#  Plot the mesh.
#
  label = ( 'step7 mesh (%d)' % ( my_grid ) )
  plot ( my_mesh, title = label )
  filename = ( 'step7_mesh_%d.png' % ( my_grid ) )
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Define the diffusivity function kappa(x).
#
  kappa = Constant ( 1.0 )
#
#  Define expressions for the exact solution and right hand side.
#
  if ( my_case == 1 ):
    u_exact_expr = Expression ( "x[0]*(1-x[0])*x[1]*(1-x[1])", degree = 5 )
    f_expr = Expression ( "2*x[1]*(1-x[1])+2*x[0]*(1-x[0])", degree = 5 )
  elif ( my_case == 2 ):
    u_exact_expr = Expression ( "exp(pow(x[0],2)+pow(x[1],2))*sin(x[0]+x[1])", degree = 10 )
    f_expr = Expression ( "-2*exp(pow(x[0],2)+pow(x[1],2))*(sin(x[0]+x[1])*(1+2*(pow(x[0],2)+pow(x[1],2))) + 2*cos(x[0]+x[1])*(x[0]+x[1]))", degree = 10 )
  elif ( my_case == 3 ):
    u_exact_expr = Expression ( "2*(1+x[1])/(pow(3+x[0],2)+pow(1+x[1],2))", degree = 10 )
    f_expr = Expression ( "0", degree = 0 )
  else:
    raise Exception ( '  The value of my_case should be between 1 and 3, your value was: {}'.format ( my_case ) )
#
#  Set spaces:
#    Es = Error estimator, 
#    Us = primal variable, 
#    Qs = interfacial flux.
#
  Es = FunctionSpace ( my_mesh, "DG",  my_degree + 2 )
  Us = FunctionSpace ( my_mesh, "CG",  my_degree + 1 )
  Qs = FunctionSpace ( my_mesh, "BDM", 1, restriction = "facet" )
#
#  Set elements, needed for the MixedElement command:
#    Ee = Error estimator element, 
#    Ue = primal variable element, 
#    Qe = interfacial flux element.
#
  Ee = FiniteElement ( "DG",  triangle, my_degree + 2 )
  Ue = FiniteElement ( "CG",  triangle, my_degree + 1 )
  Qe = FiniteElement ( "BDM", triangle, 1 )
#
#  Define the mixed element, and the corresponding function space.
#
  EUQe = MixedElement ( Ee, Ue, Qe )
  EUQs = FunctionSpace ( my_mesh, EUQe )
#
#  Extract the individual trial and test function factors from the space.
#  Here, "e" for example is a symbol for typical trial functions from Es.
#
  ( e,  u,  q  ) = TrialFunctions ( EUQs )
  ( et, ut, qt ) = TestFunctions ( EUQs )
#
#  Compute the normal vectors.
#
  n = FacetNormal ( my_mesh )
#
#  Set the components of the saddle point problem:
#
#  Define the inner product for the error estimator space.
#
  yip = dot ( grad ( e ), grad ( et ) ) * dx + e * et * dx
#
#  Set up the saddle point problem:
#
#  b( (u,q), et ) = ( kappa * grad u, grad et ) - < q.n, et >
#  This is an equation for u and q.
#
  b1 = dot ( kappa * grad ( u ), grad ( et ) ) * dx \
     - dot ( q ( '+' ), n ( '+' ) ) * ( et ( '+' ) - et ( '-' ) ) * dS \
     - dot ( q, n ) * et * ds
#
#  b( (ut,qt), e ) = ( kappa * grad ut, grad e ) - < qt.n, e >
#  This is an equation for e.
#
  b2 = dot ( kappa * grad ( ut ), grad ( e ) ) * dx \
     - dot ( qt ( '+' ), n ( '+' ) ) * ( e ( '+' ) - e ( '-' ) ) * dS \
     - dot ( qt, n ) * e * ds
#
#  Form the saddle point problem:
#
#    yip + b1 = f * et * dx
#    b2       = 0
#
  a = yip + b1 + b2 
  b = f_expr * et * dx
#
#  Indicate the Dirichlet boundary condition for the second component (u).
#
  bc = DirichletBC ( EUQs.sub(1), u_exact_expr, DomainBoundary() )
#
#  Solve the saddle point problem.
#
  euq = Function ( EUQs ) 
  solve ( a == b, euq, bcs = bc )
#
#  Extract the three components of the solution.
#  This "e, u, q" is DIFFERENT from the e, u, q above!
#  Here, "e" is the actual finite element solution function.
#
  e, u, q = euq.split ( deepcopy = True )
#
#  Plot the solution u.
#
  label = ( 'step7 u (%d)' % ( my_grid ) )
  fig = plot ( u, title = label )
  plt.colorbar ( fig )
  filename = ( 'step7_u_%d.png' % ( my_grid ) )
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Create a mesh function from the expression for the exact solution.
#
  u_exact = interpolate ( u_exact_expr, Us )
#
#  Plot the exact solution u_exact.
#
  label = ( 'step7 u_exact (%d)' % ( my_grid ) )
  fig = plot ( u_exact, title = label )
  plt.colorbar ( fig )
  filename = ( 'step7_u_exact_%d.png' % ( my_grid ) )
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Plot the error indicators.
#
  label = ( 'step7 indicators (%d)' % ( my_grid ) )
  fig = plot ( e, title = label )
  plt.colorbar ( fig )
  filename = ( 'step7_indicators_%d.png' % ( my_grid ) )
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Compute a refined mesh.
#
  my_mesh_refined = refine ( refine ( my_mesh ) )
#
#  Compute the error norm.
#
#  Warning: errornorm is NOT "symmetric" with respect to the first two arguments!
#  The computed solution must be the second argument.
#  The first argument can be an expression or a mesh function.
#
  er = errornorm ( u_exact_expr, u, norm_type = 'H1', degree_rise = 3, mesh = my_mesh_refined )
  print ( '  H1 norm of u - u_exact = %15.12f' % ( er ) )
#
#  Compute the error estimator norm.
#
  ee = sqrt ( assemble ( ( dot ( grad ( e ), grad ( e ) ) + e * e ) * dx ) )
  print ( '  Error estimator norm   = %15.12f' % ( ee ) )
#
#  Compute the error norm / error indicator ratio.
#
  print ( '  Ratio                  = %15.12f' % ( ee / er ) )
#
#  Terminate.
#
  return

def step7_test ( ):

#*****************************************************************************80
#
## step7_test tests step7.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    31 October 2018
#
#  Author:
#
#    John Burkardt
#
  import time

  print ( time.ctime ( time.time() ) )
  print ( '' )
  print ( 'step7_test:' )
  print ( '  FENICS/Python version' )
  print ( '  p-Laplacian investigation.' )
  print ( '  Exact solution known.' )
  print ( '  Compare ||u-uh|| and DPG error indicator norm.' )
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )
#
  my_case = 3

  for my_grid in [ 4, 8, 16 ]:
    my_degree = 1
    step7 ( my_grid, my_degree, my_case )
#
#  Terminate.
#
  print ( '' )
  print ( 'step7_test:' );
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )

if ( __name__ == '__main__' ):

  step7_test ( )
