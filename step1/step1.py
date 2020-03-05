#! /usr/bin/env python
#
from fenics import *

def step1 ( my_grid, my_degree ):

#*****************************************************************************80
#
## step1 solves a Poisson equation with diffusivity K(X,Y).
#
#  Discussion:
#
#    STEP1:
#      Use standard finite element approach to solve a Poisson equation
#      with a diffusivity function K(X).
#
#      - div k grad u = f in Omega = [0,1]x[0,1]
#                  u = 0 on dOmega
#
#      Use:
#        u(x,y) = sin(pi*x) * sin(pi*y)
#        k(x,y) = 1
#        f(x,y) = 2*pi*pi*sin(pi*x) * sin(pi*y)
#
#  Modified:
#
#    23 October 2018
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
#    Fluxes:           MY_DEGREE
#
  import matplotlib.pyplot as plt
#
#  Report input.
#
  print ( '' )
  print ( '  Case of %d x %d mesh' % ( my_grid, my_grid ) )
  print ( '  Approximation degree %d:' % ( my_degree ) )
#
#  Mesh the unit square.
#
  my_mesh = UnitSquareMesh ( my_grid, my_grid )
#
#  Plot the mesh.
#
  plot ( my_mesh, title = 'step1 Mesh' )
  filename = 'step1_mesh.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Define the exact solution of the PDE, U(X,Y):
#  * from which we can determine the right hand side;
#  * for convenience in setting boundary conditions;
#  * to allow for error calculations.
#
  U = Expression ( 'sin ( pi * x[0] ) * sin ( pi * x[1] )', degree = 10 )
#
#  Define the diffusivity function K(X,Y):
#
  K = Expression ( "1.0", degree = 10 )
#
#  Define the right hand side F(X,Y):
#
  F = Expression ( '2 * pi * pi * sin ( pi * x[0] ) * sin ( pi * x[1] )', degree = 10 )
#
#  Set function spaces for Error estimator, Primal variable, interfacial flux.
#
  Es = FunctionSpace ( my_mesh, 'DG',  my_degree + 2 )
  Us = FunctionSpace ( my_mesh, 'CG',  my_degree + 1 )
  Qs = FunctionSpace ( my_mesh, 'BDM', my_degree, restriction = 'facet' )
#
#  Set elements for Error estimator, Primal variable, interfacial flux.
#
  Ee = FiniteElement ( 'DG',  triangle, my_degree + 2 )
  Ue = FiniteElement ( 'CG',  triangle, my_degree + 1 )
  Qe = FiniteElement ( 'BDM', triangle, my_degree )
#
#  Define the mixed element, and the corresponding function space.
#
  EUQe = MixedElement ( Ee, Ue, Qe )
  EUQs = FunctionSpace ( my_mesh, EUQe )
#
#  Extract the individual trial and test function factors from the space.
#  Here, "e" for example is a symbol for typical trial functions from Es.
#
  ( e, u, q ) = TrialFunctions ( EUQs )
  ( et, ut, qt ) = TestFunctions ( EUQs )
#
#  Compute the normal vectors.
#
  n = FacetNormal ( my_mesh )
#
#  Define the inner product for the error estimator space.
#
  yip = dot ( grad ( e ), grad ( et ) ) * dx + e * et * dx
#
#  Set up the saddle point problem:
#
#  b( (u,q), et ) = ( grad u, grad et ) - < q.n, et >
#  This is an equation for U and Q.
#
  b1 = dot ( grad ( u ), grad ( et ) ) * dx \
     - dot ( q ( '+' ), n ( '+' ) ) * ( et ( '+' ) - et ( '-' ) ) * dS \
     - dot ( q, n ) * et * ds
#
#  b( (ut,qt), e ) = ( grad ut, grad e ) - < qt.n, e >
#  This is an equation for E.
#
  b2 = dot ( grad ( ut ), grad ( e ) ) * dx \
     - dot ( qt ( '+' ), n ( '+' ) ) * ( e ( '+' ) - e ( '-' ) ) * dS \
     - dot ( qt, n ) * e * ds
#
#  Mixed form of DPG and rhs.
#  This is an equation "indexed" by et, ut, qt.
#
  a = yip + b1 + b2 
  b = F * et * dx
#
#  Apply the Dirichlet boundary condition to the second component (U).
#
  bc = DirichletBC ( EUQs.sub(1), U, DomainBoundary() )
#
#  Solve.
#  This "e, u, q" is DIFFERENT from the e, u, q above!
#  Here, "e" is the actual finite element solution function.
#
  euq = Function ( EUQs ) 
  solve ( a == b, euq, bcs = bc )
  e, u, q = euq.split ( deepcopy = True )
#
#  Plot the solution.
#
  plot ( u, title = 'step1 solution' )
  filename = 'step1_solution.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Plot the exact solution.
#
  plot ( U, mesh = my_mesh, title = 'step1 exact' )
  filename = 'step1_exact.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Norm of error over given mesh.
#
  er = errornorm ( U, u, norm_type = 'L2', degree_rise = 3, mesh = my_mesh )
  print ( '  L2-norm of  (u - uh)   = %15.12f' % er )
#
#  Refine the mesh.
#
  my_fine_mesh = refine ( refine ( my_mesh ) )
#
#  Plot the refined mesh.
#
  plot ( my_fine_mesh, title = 'step1 Refined Mesh' )
  filename = 'step1_refined_mesh.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Compare U (exact solution expression) and u (computed solution)
#
  print ( '' )
  print ( '  Use refined mesh for further error analysis.' )
  print ( '' )
  er = errornorm ( U, u, norm_type = 'H1', degree_rise = 3, mesh = my_fine_mesh )
  print ( '  H1-norm of  (u - uh)   = %15.12f' % er )
  ee = sqrt ( assemble ( ( dot ( grad ( e ), grad ( e ) ) + e * e ) * dx ) )
  print ( '  Error estimator norm   = %15.12f' % ee )
#
#  Compute the H1 projection of U (a standard Galerkin solution)
#
  uu = TrialFunction ( Us )
  vv = TestFunction ( Us )
  aa = ( dot ( grad ( uu ), grad ( vv ) ) + uu * vv ) * dx
  bb = ( F + U ) * vv * dx
  bc = DirichletBC ( Us, U, DomainBoundary() )
  up = Function ( Us )
  solve ( aa == bb, up, bcs = bc )
#
#  Compare U (exact solution expression) and up (H1 projection)
#
  erp = errornorm ( U, up, norm_type = 'H1', degree_rise = 3, mesh = my_fine_mesh )
  print ( '  Error in H1 Projection = %15.12f' % ( erp ) )
  if ( 1.0e-12 < abs ( erp ) ):
    print ( '  Error:Projection ratio = %15.12f' % ( er / erp ) )
  print ( '  H1-norm of (uh - up) = %15.12f' % \
    errornorm ( u, up, norm_type = 'H1', degree_rise = 3 ) )
#
#  Terminate.
#
  return

def step1_test ( ):

#*****************************************************************************80
#
## step1_test tests step1.
#
#  Modified:
#
#    23 October 2018
#
#  Author:
#
#    John Burkardt
#
  import time

  print ( time.ctime ( time.time() ) )
  print ( '' )
  print ( 'step1:' )
  print ( '  FENICS/Python version' )
  print ( '  First step in P-Laplacian investigation.' )
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  my_grid = 4
  my_degree = 1
  step1 ( my_grid, my_degree )
#
#  Terminate.
#
  print ( '' )
  print ( 'step1:' );
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )

if ( __name__ == '__main__' ):

  step1_test ( )
