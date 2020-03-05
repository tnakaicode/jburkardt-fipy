#! /usr/bin/env python
#
from fenics import *

def p_laplacian ( my_grid, my_degree, my_p ):

#*****************************************************************************80
#
## p_laplacian tries to solve the p-Laplacian problem directly.
#
#  Discussion:
#
#    Solve the nonlinear PDE:
#
#      - div gamma(u) grad u = f in Omega = [0,1]x[0,1]
#                          u = 0 on dOmega
#
#    where:
#
#      gamma(u) = ( epsilon^2 + |grad u|^2) ) ^ ( (p-2)/2 )
#               = ( epsilon^2 + wx^2 + wy^2 ) ^ ( (p-2)/2 )
#
#      epsilon = 1.0E-05
#
#      f(x,y) = 2*pi*pi*sin(pi*x) * sin(pi*y)
#
#    The parameter epsilon artificially keeps gamma(u) bounded above 0,
#    guaranteeing that the power (p-2)/2 is meaningful, and the equation
#    is solvable.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 October 2018
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    John Barrett, W B Liu,
#    Finite Element Approximation of the p-Laplacian,
#    Mathematics of Computation,
#    Volume 61, Number 204, pages 423-537, October 1993.
#
#    Patrick Farrell, Hans Petter Langtangen, Marie Rognes, Garth Wells,
#    Implementing finite element models in Python/FEniCS:
#    steady nonlinear PDEs,
#    MMSC: Python in Scientific Computing,
#    May 6, 2015
#
#  Parameters:
#
#    Input, integer MY_GRID, the number of (pairs of) elements to use in the 
#    X and Y directions.  Number of triangular elements is 2*my_grid*my_grid.
#
#    Input, integer MY_DEGREE, specifies approximation degree for the solution.
#
#    Input, real MY_P, the desired value of the p exponent in the equation.
#    2 represents the usual Laplacian.
#    0 < P < 2 will involve dividing by the norm of the gradient, which will
#    result in a large coefficient anywhere near a solution extremum,
#    and (if epsilon were 0) an unbounded value at an extremum.
#
  import matplotlib.pyplot as plt
  import numpy as np
#
#  Report input.
#
  print ( '' )
  print ( '  Case of %d x %d mesh' % ( my_grid, my_grid ) )
  print ( '  Approximation degree %d:' % ( my_degree ) )
  print ( '  P-Laplacian parameter P = %g' % ( my_p ) )
#
#  Mesh the unit square.
#  If we don't use the "crossed" option, then two corner triangles are entirely 0.
#
  my_mesh = UnitSquareMesh ( my_grid, my_grid, 'crossed' )
#
#  Plot the mesh.
#
  plot ( my_mesh, title = 'p_laplacian Mesh' )
  filename = 'p_laplacian_mesh.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Define the function space.
#
  V = FunctionSpace ( my_mesh, "Lagrange", my_degree )
#
#  Because this is a nonlinear problem, u is a function, not a trial function!
#  We need to give it an initial value for the Newton iteration.
#  The function here is positive in the interior, and satisfies the boundary conditions.
#
  uexpr = Expression ( "x[0] * ( 1.0 - x[0] ) * x[1] * ( 1.0 - x[1] )", degree = 10 )
  u = interpolate ( uexpr, V )
#
#  v is just a test function.
#
  v = TestFunction ( V )
#
#  Define the right hand side f(x,y):
#
# f = Expression ( '2 * pi * pi * sin ( pi * x[0] ) * sin ( pi * x[1] )', degree = 10 )
  f = Constant ( 1.0 )
#
#  Define the boundary condition.
#
  g = Constant ( 0.0 )
  bc = DirichletBC ( V, g, DomainBoundary ( ) )
#
#  Define the "diffusivity".
#
  epsilon = Constant ( 1.0e-05 )
  p = Constant ( my_p )

  def gamma ( u ):
    value = ( epsilon ** 2 + inner ( grad ( u ), grad ( u ) ) ) ** ( ( p - 2 ) / 2 )
    return value
#
#  Define the nonlinear function.
#
  F = inner ( grad ( v ), gamma ( u ) * grad ( u ) ) * dx - inner ( f, v ) * dx
#
#  Request a solution.
#
  solve ( F == 0, u, bc )
#
#  Plot the solution u.
#
  fig = plot ( u, title = 'p_laplacian solution' )
  plt.colorbar ( fig )
  filename = 'p_laplacian_solution.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Plot |grad u|.
#
  gradun = sqrt ( grad ( u ) ** 2 )
  fig = plot ( gradun, title = 'p_laplacian |grad u|' )
  plt.colorbar ( fig )
  filename = 'p_laplacian_gradun.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Plot gamma(u).
#
  gammau = gamma ( u )
  fig = plot ( gammau, title = 'p_laplacian gamma(u)' )
  plt.colorbar ( fig )
  filename = 'p_laplacian_gammau.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Terminate.
#
  return

def p_laplacian_test ( ):

#*****************************************************************************80
#
## p_laplacian_test tests p_laplacian.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    29 October 2018
#
#  Author:
#
#    John Burkardt
#
  import time

  print ( time.ctime ( time.time() ) )
  print ( '' )
  print ( 'p_laplacian:' )
  print ( '  FENICS/Python version' )
  print ( '  Step 6 in p-Laplacian investigation.' )
  print ( '  Solve the nonlinear equation directly.' )
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )
#
#  Set input parameters:
#  P = 3.0 means (P-2)/2 is positive.
#  P = 1.5 means (P-2)/2 is negative.
#
  my_grid = 16
  my_degree = 1
  my_p = 1.5

  p_laplacian ( my_grid, my_degree, my_p )
#
#  Terminate.
#
  print ( '' )
  print ( 'p_laplacian:' );
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )

if ( __name__ == '__main__' ):

  p_laplacian_test ( )
