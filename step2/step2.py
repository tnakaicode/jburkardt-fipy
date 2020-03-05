#! /usr/bin/env python
#
from fenics import *

def step2 ( my_grid, my_degree ):

#*****************************************************************************80
#
## step2 solves a Poisson equation with piecewise constant diffusivity K(X,Y).
#
#  Discussion:
#
#    Use the mixed DPG method to solve a Poisson equation
#    with a piecewise constant diffusivity function K(X).
#
#      - div k grad u = f in Omega = [0,1]x[0,1]
#                  u = 0 on dOmega
#
#    Use:
#      k(x,y) = 1 if x < 1/3
#             = 4 otherwise
#      f(x,y) = 2*pi*pi*sin(pi*x) * sin(pi*y)
#
#    For this case, I could have used an "Expression()" for k(x,y)
#    involving the C ternary operator:
#      k(x,y) = condition ? true value : false value
#    but this can't be easily generalized to more complicated cases.
#
#    The step1 problem used k(x,y)=1.  Since step2 uses a changed k,  
#    the exact solution u is not known.  However, we retain the
#    source term f(x,y) and zero boundary conditions, for lack of any 
#    better case to consider.
#
#  Modified:
#
#    24 October 2018
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
  plot ( my_mesh, title = 'step2 Mesh' )
  filename = 'step2_mesh.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Define the diffusivity function K(X,Y).
#  The complication here is that the input "x" is not a Python numerical
#  array, but a corresponding, but symbolic, object.  So we have to use
#  the UFL "conditional" function to express a test on X.
#
  def diffusivity ( x ):
#   value = conditional ( le ( 3.0 * x[0], 1.0 ), 1.0, 4.0 )
    value = conditional ( le ( x[0] + x[1], 1.0 ), 1.0, 4.0 )
    return value
#
#  Plot the diffusivity function.
#
  x = SpatialCoordinate ( my_mesh )

  plot ( diffusivity ( x ), title = 'Step2 Diffusivity' )
  filename = 'step2_diffusivity.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
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
#  b( (u,q), et ) = ( k * grad u, grad et ) - < q.n, et >
#  This is an equation for U and Q.
#
  b1 = dot ( diffusivity ( x ) * grad ( u ), grad ( et ) ) * dx \
     - dot ( q ( '+' ), n ( '+' ) ) * ( et ( '+' ) - et ( '-' ) ) * dS \
     - dot ( q, n ) * et * ds
#
#  b( (ut,qt), e ) = ( k * grad ut, grad e ) - < qt.n, e >
#  This is an equation for E.
#
  b2 = dot ( diffusivity ( x ) * grad ( ut ), grad ( e ) ) * dx \
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
#
  euq = Function ( EUQs ) 
  solve ( a == b, euq, bcs = bc )
  e, u, q = euq.split ( deepcopy = True )
#
#  Plot the solution.
#
  fig = plot ( u, title = 'step2 solution' )
  plt.colorbar ( fig )
  filename = 'step2_solution.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Plot the error indicators.
#
  fig = plot ( e, title = 'step2 indicators' )
  plt.colorbar ( fig )
  filename = 'step2_indicators.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Terminate.
#
  return

def step2_test ( ):

#*****************************************************************************80
#
## step2_test tests step2.
#
#  Modified:
#
#    24 October 2018
#
#  Author:
#
#    John Burkardt
#
  import time

  print ( time.ctime ( time.time() ) )
  print ( '' )
  print ( 'step2:' )
  print ( '  FENICS/Python version' )
  print ( '  Step 2 in P-Laplacian investigation.' )
  print ( '  Piecewise constant diffusivity.' )
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  my_grid = 8
  my_degree = 1
  step2 ( my_grid, my_degree )
#
#  Terminate.
#
  print ( '' )
  print ( 'step2:' );
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )

if ( __name__ == '__main__' ):

  step2_test ( )
