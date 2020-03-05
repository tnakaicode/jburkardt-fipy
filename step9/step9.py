#! /usr/bin/env python
#
from fenics import *

def step9 ( my_grid, my_degree, my_case ):

#*****************************************************************************80
#
## step9 plots error behavior for adapative refinement cases.
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
#    Do this for a sequence of refined meshes, and at the end,
#    plot the number of elements and the estimated errors.
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
#    if my_case == 4:
#
#      u_exact(x,y) = exp ( -x^2-y^2 )
#      f(x,y)       = 4 * ( - 1 + x^2 + y^2 ) * exp ( -x^2-y^2 )
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
#    1 <= MY_CASE <= 4.
#
  import matplotlib.pyplot as plt
  import numpy as np
#
#  Report input.
#
  REFINE_RATIO = 0.50

  print ( '' )
  print ( '  Exact case #%d' % ( my_case ) )
  print ( '  %dx%d mesh on the unit square' % ( my_grid, my_grid ) )
  print ( '  Approximation degree %d:' % ( my_degree ) )
  print ( '  Proportion of cells to refine is %g' % ( REFINE_RATIO ) )
#
#  Internal parameters:
#    REFINE_RATIO: proportion of cells to refine.
#    LEVEL_MAX: maximum number of iterations
#
  LEVEL_MAX = 5
#
#  Mesh the unit square.
#
  my_mesh = UnitSquareMesh ( my_grid, my_grid, 'right/left' )
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
  elif ( my_case == 4 ):
    u_exact_expr = Expression ( "exp(-pow(x[0],2)-pow(x[1],2))", degree = 10 )
    f_expr = Expression ( "-(-4*pow(x[0],2)-4*pow(x[1],2)+4)*exp(-pow(x[0],2)-pow(x[1],2))", degree = 10 )
  else:
    raise Exception ( '  The value of my_case should be between 1 and 4, your value was: {}'.format ( my_case ) )
#
#  Initialize Numpy arrays that will hold, for each level, the number of
#  cells (elements) and the error estimate.
#
  nvec = np.zeros ( 0 )
  evec = np.zeros ( 0 )
#
#  ADAPTIVE LOOP
#
  for level in range ( LEVEL_MAX ):
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
#  Compute the error estimate.
#
    cell_num = my_mesh.num_cells ( )
    error_estimate = sqrt ( assemble ( ( dot ( grad ( e ), grad ( e ) ) + e * e ) * dx ) )
    print ( '  Level %d: %d cells, error estimate is %g' % ( level, cell_num, error_estimate ) )
#
#  Add the current number of cells and error estimate to our lists.
#
    nvec = np.append ( nvec, cell_num )
    evec = np.append ( evec, error_estimate )
#
#  If this is the last level, we're done now.
#  Otherwise, compute error indicators and refine before next level.
#
    if ( LEVEL_MAX - 1 <= level ):
      print ( '' )
      print ( '  Finished with adaptive refinement sequence.' )
      break
#
#  Hs: an piecewise constant function space for error indicators.
#
    hs_degree = 0
    Hs = FunctionSpace ( my_mesh, 'DG', hs_degree )
#
#  Define a test function.
#
    ht = TestFunction ( Hs )
#
#  Define g, which is a piecewise constant mesh function whose value in each
#  cell is ||e||^2
#
    g = assemble ( ( dot ( grad ( e ), grad ( e ) ) + e * e ) * ht * dx )
#
#  Replace the mesh function by an array.
#
    g = g.get_local ( )
#
#  Create a cell marker array.
#
    cell_markers = MeshFunction ( 'bool', my_mesh, my_mesh.topology().dim() )
#
#  This strange ugly command descending sorts G, then indexes the value that is
#  halfway down, and calls that value g0.  
#
    g0_index = int ( len ( g ) * REFINE_RATIO )
    print ( '  Refining %d cells' % ( g0_index + 1 ) )

    g0 = sorted ( g, reverse = True ) [ int ( len ( g ) * REFINE_RATIO ) ]
#
#  Mark cells with a high error norm.
#
    for c in cells ( my_mesh ):
      cell_markers[c] = g[c.index()] > g0
#
#  Refine the mesh.
#
    my_mesh = refine ( my_mesh, cell_markers )
#
#  Report N, Error Estimate
#
  print ( '' )
  print ( '  N  Error estimate:' )
  print ( '' )
  for level in range ( 0, LEVEL_MAX ):
    print ( '  %6d  %10.4g' % ( nvec[level], evec[level] ) )
#
#  Convert our data to logarithms base 10.
#
  nvec = np.log10 ( nvec )
  evec = np.log10 ( evec )
#
#  Plot error behavior.
#
  plt.plot ( nvec, evec, marker = "o" )
  plt.title ( 'step9 error decay case %d' % ( my_case ) )
  plt.grid ( True )
  plt.xlabel ( '<-- Log10 ( Cells ) -->' )
  plt.ylabel ( '<-- Log10 ( Error ) -->' )
  filename = ( 'step9_error_decay_case%d.png' % ( my_case ) )
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Terminate.
#
  return

def step9_test ( ):

#*****************************************************************************80
#
## step9_test tests step9.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 November 2018
#
#  Author:
#
#    John Burkardt
#
  import time

  print ( time.ctime ( time.time() ) )
  print ( '' )
  print ( 'step9_test:' )
  print ( '  FENICS/Python version' )
  print ( '  p-Laplacian investigation.' )
  print ( '  Display error decay for adaptive refinement.' )
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  for my_case in range ( 1, 5 ):
    my_grid = 4
    my_degree = 2
    step9 ( my_grid, my_degree, my_case )
#
#  Terminate.
#
  print ( '' )
  print ( 'step9_test:' );
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )

if ( __name__ == '__main__' ):

  step9_test ( )
