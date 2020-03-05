#! /usr/bin/env python
#
from fenics import *

def step5 ( my_grid, my_degree ):

#*****************************************************************************80
#
## step5 solves a Poisson equation with diffusivity |grad(W)|.
#
#  Discussion:
#
#    Use the mixed DPG method to solve a Poisson equation
#    with a piecewise constant diffusivity function K(X).
#
#      - div k(x,y) grad u = f in Omega = [0,1]x[0,1]
#                        u = 0 on dOmega
#
#    Use:
#      k(x,y) = |grad w| = sqrt ( wx^2 + wy^2 )
#      f(x,y) = 2*pi*pi*sin(pi*x) * sin(pi*y)
#    Where:
#      w(x,y) is, for right now, an arbitrary scalar field, which we 
#      will construct as the solution of
#      - div grad w = f in Omega = [0,1]x[0,1]
#                 w = 0 on dOmega
#
#    Later, we will be working with an iterative procedure, starting with
#    u0, and producing u1, u2, ... and w(x,y) will be the value of u
#    from the previous iterative step.
#
#  Modified:
#
#    26 October 2018
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
  import numpy as np
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
  plot ( my_mesh, title = 'step5 Mesh' )
  filename = 'step5_mesh.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#-------------------------------------------------------------------------------
#  Set up and solve for a field W.
#  |grad(w)| will be used as the diffusivity for U.
#-------------------------------------------------------------------------------
  F = Expression ( '2 * pi * pi * sin ( pi * x[0] ) * sin ( pi * x[1] )', degree = 10 )
  Ws = FunctionSpace ( my_mesh, 'CG',  my_degree + 1 )
  wt = TrialFunction ( Ws )
  vt = TestFunction ( Ws )
  awv = dot ( grad ( wt ), grad ( vt ) ) * dx
  fv = F * vt * dx
  w = Function ( Ws )
  wbc = DirichletBC ( Ws, w, DomainBoundary() )
  solve ( awv == fv, w, wbc )
#
#  Plot w.
#
  plot ( w, title = 'step5 Diffusivity' )
  filename = 'step5_diffusivity.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )

#
#  Plot w.
#
  plot ( grad ( w )**2, title = 'step5 Diffusivity grad norm' )
  filename = 'step5_diffusivity_grad_norm.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )

#-------------------------------------------------------------------------------
#  Now we seek to solve
#    - div |grad(W)| grad U = F
#-------------------------------------------------------------------------------
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
#  b( (u,q), et ) = ( k * grad u, grad et ) - < q.n, et >
#  This is an equation for U and Q.
#
  b1 = dot ( sqrt ( grad ( w )**2 ) * grad ( u ), grad ( et ) ) * dx \
     - dot ( q ( '+' ), n ( '+' ) ) * ( et ( '+' ) - et ( '-' ) ) * dS \
     - dot ( q, n ) * et * ds
#
#  b( (ut,qt), e ) = ( k * grad ut, grad e ) - < qt.n, e >
#  This is an equation for E.
#
  b2 = dot ( sqrt ( grad ( w )**2 ) * grad ( ut ), grad ( e ) ) * dx \
     - dot ( qt ( '+' ), n ( '+' ) ) * ( e ( '+' ) - e ( '-' ) ) * dS \
     - dot ( qt, n ) * e * ds
#
#  Set the saddle point problem:
#
#    yip + b1 = F * et * dx
#    b2       = 0
#
  a = yip + b1 + b2
  b = F * et * dx
#
#  Apply the Dirichlet boundary condition to the second component (U).
#
  U = Expression ( "0", degree = 10 )
  bc = DirichletBC ( EUQs.sub(1), U, DomainBoundary() )
#
#  Solve.
#  This "e, u, q" is DIFFERENT from the e, u, q above!
#  Here, "e" is the actual finite element solution function.
#  Specify a nonlinear equation to be solved.
#
  euq = Function ( EUQs ) 
  solve ( a == b, euq, bcs = bc )
  e, u, q = euq.split ( deepcopy = True )
#
#  Plot the solution u.
#
  fig = plot ( u, title = 'step5 solution' )
  plt.colorbar ( fig )
  filename = 'step5_solution.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Plot |grad u|.
#
  gradun = sqrt ( grad ( u ) ** 2 )
  fig = plot ( gradun, title = 'step5 |grad u|' )
  plt.colorbar ( fig )
  filename = 'step5_gradun.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Plot the error indicators.
#
  fig = plot ( e, title = 'step5 indicators' )
  plt.colorbar ( fig )
  filename = 'step5_indicators.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Terminate.
#
  return

def step5_test ( ):

#*****************************************************************************80
#
## step5_test tests step5.
#
#  Modified:
#
#    25 October 2018
#
#  Author:
#
#    John Burkardt
#
  import time

  print ( time.ctime ( time.time() ) )
  print ( '' )
  print ( 'step5:' )
  print ( '  FENICS/Python version' )
  print ( '  Step 5 in P-Laplacian investigation.' )
  print ( '  Diffusivity is |grad(u)|.' )
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  my_grid = 16
  my_degree = 1
  step5 ( my_grid, my_degree )
#
#  Terminate.
#
  print ( '' )
  print ( 'step5:' );
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )

if ( __name__ == '__main__' ):

  step5_test ( )
