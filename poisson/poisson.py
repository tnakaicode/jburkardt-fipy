#! /usr/bin/env python
#
from fenics import *

def poisson ( ):

#*****************************************************************************80
#
## poisson solves the Poisson problem on the unit square.
#
#  Discussion:
#
#    Del^2 U = - 6             in Omega
#    U(dOmega) = Uexact(x,y)   on dOmega
#
#    Uexact = 1 + x^2 + 2 y^2
#    Omega = unit square [0,1]x[0,1]
#
#  Modified:
#
#    21 October 2018
#
#  Author:
#
#    John Burkardt
#
#  Reference:
#
#    Hans Petter Langtangen, Anders Logg,
#    Solving PDEs in Python - The FEniCS Tutorial Volume !
#
  import matplotlib.pyplot as plt
#
#  Create an 8x8 triangular mesh on the unit square.
#
  mesh = UnitSquareMesh ( 8, 8 )
#
#  Define the function space.
#
  V = FunctionSpace ( mesh, 'P', 1 )
#
#  Define the exact solution.
#  Setting "degree = 2" should be OK.
#  FENICS prefers degree = 4 for accuracy...
#
  u_D = Expression ( '1 + x[0]*x[0] + 6*x[1]*x[1]', degree = 4 )
#
#  Define the boundary condition.
#
  def boundary ( x, on_boundary ):
    return on_boundary

  bc = DirichletBC ( V, u_D, boundary )
#
#  Define the variational problem.
#
  u = TrialFunction ( V )
  v = TestFunction ( V )
  f = Constant ( -14.0 )
  a = dot ( grad ( u ), grad ( v ) ) * dx
  L = f * v * dx
#
#  Compute the solution.
#
  u = Function ( V )
  solve ( a == L, u, bc )
#
#  Plot the mesh.
#
  plot ( mesh, title = 'Mesh for Poisson equation' )
  filename = 'poisson_mesh.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Plot the solution.
#
  plot ( u, mode = 'contour', title = 'Solution for Poisson equation' )
  filename = 'poisson_solution.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Compute the error in the L2 norm.
#
  error_L2 = errornorm ( u_D, u, 'L2' )
  print ( 'error_L2  =', error_L2 )
#
#  Compute maximum error at vertices.
#
  vertex_values_u_D = u_D.compute_vertex_values ( mesh )
  vertex_values_u = u.compute_vertex_values ( mesh )
  import numpy as np
  error_max = np.max ( np.abs ( vertex_values_u_D - vertex_values_u ) )
  print ( 'error_max =', error_max )

  return

def poisson_test ( ):

#*****************************************************************************80
#
## poisson_test tests poisson.
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
  print ( 'poisson_test:' )
  print ( '  FENICS/Python version' )
  print ( '  Poisson equation on the unit square.' )

  poisson ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'ell_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  poisson_test ( )

