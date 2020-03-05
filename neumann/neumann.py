#! /usr/bin/env python
#
from fenics import *

def neumann ( ):

#*****************************************************************************80
#
## neumann solves a boundary value problem with homogenous Neumann conditions.
#
#  Discussion:
#
#    - div grad u + u = f        in Omega
#    du/dn(dOmega) = 0           on dOmega
#
#    u_exact = (x*(x-2)*sin(2*pi*y))^2
#    f = -8 * (pi*x*(x-2))^2 * cos(4*pi*y)
#        + (x^4-4x^3-8x^2+24x-8) * (sin(2*pi*y))^2
#    Omega = unit square [0,1]x[0,1]
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    07 November 2018
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Doug Arnold
#    Getting Started with Fenics
#
  import matplotlib.pyplot as plt
#
#  Create an 8x8 triangular mesh on the unit square.
#
  mesh = UnitSquareMesh ( 10, 20 )
#
#  Define the function space.
#
  Vh = FunctionSpace ( mesh, 'Lagrange', 2 )
#
#  Define the exact solution.
#
  u_exact_expr = Expression ( 'pow(x[0]*(x[0]-2.0)*sin(2*pi*x[1]),2)', degree = 10 )
#
#  Define the right hand side.
#
  f_expr = Expression ( '-8*pow(pi*(pow(x[0],2)-2*x[0]),2)*cos(4*pi*x[1]) \
    +(pow(x[0],4)-4*pow(x[0],3)-8*pow(x[0],2)+24*x[0]-8) \
    * pow(sin(2*pi*x[1]),2)', degree = 10 )
# f_expr = Constant ( 1.0 )
#
#  Define the variational problem.
#
  u = TrialFunction ( Vh )
  v = TestFunction ( Vh )
  b = ( dot ( grad ( u ), grad ( v ) ) + u * v ) * dx
  F = f_expr * v * dx
#
#  Compute the solution.
#  Because this is a homogeneous Neumann problem, no boundary conditions
#  are specified.
#
  uh = Function ( Vh )
  solve ( b == F, uh )
#
#  Plot the mesh.
#
  plot ( mesh, title = 'Neumann mesh' )
  filename = 'neumann_mesh.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Plot the solution.
#
  fig = plot ( uh, title = 'Neumann solution' )
  plt.colorbar ( fig )
  filename = 'neumann_solution.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Plot the solution.
#
  plot ( uh, mode = 'contour', title = 'Neumann contour' )
  plt.colorbar ( fig )
  filename = 'neumann_contour.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Compute the error in the L2 norm.
#
  error_L2 = errornorm ( u_exact_expr, uh, 'L2' )
  print ( '  error_L2  =', error_L2 )
#
#  Terminate.
#
  return

def neumann_test ( ):

#*****************************************************************************80
#
## neumann_test tests neumann.
#
#  Modified:
#
#    07 November 2018
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
  print ( 'neumann_test:' )
  print ( '  FENICS/Python version' )
  print ( '  Boundary value problem with Neumann conditions.' )

  neumann ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'neumann_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  neumann_test ( )
