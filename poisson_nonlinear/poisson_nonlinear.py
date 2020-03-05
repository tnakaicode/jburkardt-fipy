#! /usr/bin/env python
#
from fenics import *

def poisson_nonlinear ( ):

#*****************************************************************************80
#
## poisson_nonlinear solves the nonlinear Poisson problem on the unit square.
#
#  Discussion:
#
#    - div( ( 1 + u^2 ) grad(u) ) = x * sin(y)  in Omega
#    U = 1   if x = 1 on dOmega
#    U = 0   otherwise on dOmega
#
#    Omega = unit square [0,1]x[0,1]
#
#  Modified:
#
#    20 October 2018
#
#  Author:
#
#    John Burkardt
#
  import matplotlib.pyplot as plt
#
#  Set up the mesh.
# 
  mesh = UnitSquareMesh ( 32, 32 )
#
#  Function space.
#
  V = FunctionSpace ( mesh, "CG", 1 )
#
#  Dirichlet boundary conditions.
# 
  def right_boundary ( x, on_boundary ):
    return abs ( x[0] - 1.0 ) < DOLFIN_EPS and on_boundary

  g = Constant ( 1.0 )
  bc = DirichletBC ( V, g, right_boundary )
#
#  Set up the variational form.
# 
  u = Function ( V )
  v = TestFunction ( V )
  f = Expression ( "x[0]*sin(x[1])", degree = 10 )

  F = inner ( ( 1 + u**2 ) * grad(u), grad(v) ) * dx - f * v * dx
#
#  Compute the solution.
#  Because this is a nonlinear equation, we don't use the "auv = fv" form
#  for the solver, we ask it to solve "F==0" instead.
# 
  solve ( F == 0.0, u, bc, \
    solver_parameters = { "newton_solver":{"relative_tolerance":1e-6} } )
# 
#  Plot the solution.
# 
  plot ( u, title = "Nonlinear Poisson Solution" )
  filename = 'poisson_nonlinear_solution.png'
  plt.savefig ( filename )
  print ( 'Grpahics saved as "%s"' % ( filename ) )
  plt.close ( )
# 
#  Plot the gradient.
# 
  plot ( grad ( u ), title = "Nonlinear Poisson Gradient" )
  filename = 'poisson_nonlinear_gradient.png'
  plt.savefig ( filename )
  print ( 'Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Terminate.
#
  return

def poisson_nonlinear_test ( ):

#*****************************************************************************80
#
## poisson_nonlinear_test tests poisson_nonlinear.
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
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  print ( '' )
  print ( 'poisson_nonlinear_test:' )
  print ( '  FENICS/Python version' )
  print ( '  Solve a nonlinear Poisson problem.' )

  poisson_nonlinear ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'poisson_nonlinear_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  poisson_nonlinear_test ( )
