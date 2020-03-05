#! /usr/bin/env python3
#
from fenics import *

def heat_steady_03 ( ):

#*****************************************************************************80
#
## heat_steady_03, 2D steady heat equation on a rectangle.
#
#  Discussion:
#
#    Heat equation in a long rectangle with varying thermal diffusivity k(x,y).
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 October 2018
#
#  Author:
#
#    John Burkardt
#
  import matplotlib.pyplot as plt
#
#  Define the mesh.
#
  x_left = 0.0
  x_right = 5.0
  y_bottom = 0.0
  y_top = 5.0

  sw = Point ( x_left, y_bottom )
  ne = Point ( x_right, y_top )
  mesh = RectangleMesh ( sw, ne, 50, 50 )
#
#  Define the function space.
#
  V = FunctionSpace ( mesh, "Lagrange", 1 )
#
#  Define boundary conditions on left, and top/bottom.
#
  u_left = Expression ( "100 * x[1] * ( 5.0 - x[1] )", degree = 3 )
  def on_left ( x, on_boundary ):
    return ( x[0] <= x_left + DOLFIN_EPS )
  bc_left = DirichletBC ( V, u_left, on_left )

  u_tb = 0.0
  def on_tb ( x, on_boundary ):
    return ( x[1] <= y_bottom + DOLFIN_EPS or y_top - DOLFIN_EPS <= x[1] )
  bc_tb = DirichletBC ( V, u_tb, on_tb )

  bc = [ bc_left, bc_tb ]
#
#  Define the trial functions (u) and test functions (v).
#
  u = TrialFunction ( V )
  v = TestFunction ( V )
#
#  Write the bilinear form.
#
  k = Expression ( "1.0 + 100 * ( 5.0 - x[1] ) * x[1]", degree = 3 )
  Auv = k * inner ( grad ( u ), grad ( v ) ) * dx
#
#  Write the linear form.
#
  f = Constant ( 0.0 )
  Lv = f * v * dx
#
#  Solve the variational problem with boundary conditions.
#
  u = Function ( V )
  solve ( Auv == Lv, u, bc )
#
#  Plot the solution.
#
  plot ( u, title = 'heat_steady_03' )
  filename = 'heat_steady_03.png'
  plt.savefig ( filename )
  print ( "Saving graphics in file '%s'" % ( filename ) )
  plt.close ( )
#
#  Terminate.
#
  return

def heat_steady_03_test ( ):

#*****************************************************************************80
#
## heat_steady_03_test tests heat_steady_03.
#
#  Modified:
#
#    21 October 2018
#
#  Author:
#
#    John Burkardt
#
  import time

  print ( time.ctime ( time.time() ) )
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  print ( '' )
  print ( 'heat_steady_03_test:' )
  print ( '  FENICS/Python version' )
  print ( '  Solve the heat equation over a square' )
  print ( '  with smoothly-varying thermal diffusivity.' )

  heat_steady_03 ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'heat_steady_03_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  heat_steady_03_test ( )
