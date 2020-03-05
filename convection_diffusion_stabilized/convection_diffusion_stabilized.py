#! /usr/bin/env python
#
from fenics import *

def convection_diffusion_stabilized ( my_mu, my_grid, my_beta ):

#*****************************************************************************80
#
## convection_diffusion_stabilized simulates a 1D convection diffusion problem.
#
#  Discussion:
#
#    The equation to be solved is:
#
#      - ux - mu uxx = f in Omega = the unit interval.
#        u(0) = 0
#        u(1) = 1
#
#    where 
#
#      mu > 0
#
#      u_exact = ( e^(-x/mu) - 1 ) / ( e^(-1/mu) - 1 )
#
#      f(x) = 0
#
#    Especially as mu approaches zero, the solutions become oscillatory.
#    To stabilize the solutions, we add an artificial viscosity term,
#    controlled by a parameter beta, and multiplied by the smallest 
#    mesh spacing hmin, so that the equation becomes
#
#      - ux - (mu+beta hmin) uxx = f in Omega = the unit interval.
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
#  Reference:
#
#    Anders Logg, Kent-Andre Mardal,
#    Lectures on the Finite Element Method.
#
#  Parameters:
#
#    Input, real my_mu, the viscosity.
#    0 < my_mu.
#
#    Input, integer my_grid, the resolution on the unit interval.
#
#    Input, integer my_beta, the artificial viscosity parameter.
#    my_beta = 0 is the original problem.
#    my_beta = 0.5 is a reasonable value for stabilization.
#
  import matplotlib.pyplot as plt
#
#  Set the mesh.
#
  print ( '' )
  print ( '  Unit interval mesh n = %d' % ( my_grid ) )
  my_mesh = UnitIntervalMesh ( my_grid )
#
#  Set the function space.
#
  V = FunctionSpace ( my_mesh, 'CG', 1 )
#
#  Set the trial and test functions.
#
  u = TrialFunction ( V )
  v = TestFunction ( V )
#
#  Make a copy of my_mu.
#
  mu = Constant ( my_mu )
#
#  Set the right hand side.
#
  f = Constant ( 0.0 )
#
#  Define the minimum mesh spacing.
#
  hmin = my_mesh.hmin ( )
#
#  Make a copy of beta.
#
  beta = Constant ( my_beta )
#
#  Define the bilinear form and right hand side
#
  a = ( - u.dx(0) * v + mu * u.dx(0) * v.dx(0) + beta * hmin * u.dx(0) * v.dx(0) ) * dx
  L = f * v * dx
#
#  Define the exact solution.
#
  u_expr = Expression ( "(exp(-x[0]/%e)-1)/(exp(-1/%e)-1)" % ( my_mu, my_mu ), degree = 10 )
#
#  Define the boundary condition.
#
  def boundary ( x ):
    value = x[0] < DOLFIN_EPS or 1.0 - DOLFIN_EPS < x[0]
    return value

  bc = DirichletBC ( V, u_expr, boundary )
#
#  Solve.
#
  uh = Function ( V )
#
#  Solve the system.
#
  solve ( a == L, uh, bc )
#
#  Project the exact solution.
#
  u_exact = interpolate ( u_expr, V )
#
#  Plot the solution.
#
  fig = plt.figure ( )
  ax = plt.subplot ( 111 )
  plot (uh, label = 'Computed' )
  plot (u_exact, label = 'Exact' )
  ax.legend ( )
  ax.grid ( True )
  plt.title ( 'Convection_diffusion_stabilized solutions, grid %d' % ( my_grid ) )
  filename = ( 'convection_diffusion_stabilized_solutions_grid%d.png' % ( my_grid ) )
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Terminate.
#
  return

def convection_diffusion_stabilized_test ( ):

#*****************************************************************************80
#
## convection_diffusion_stabilized_test tests convection_diffusion_stabilized.
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
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  print ( '' )
  print ( 'convection_diffusion_stabilized_test:' )
  print ( '  FENICS/Python version' )
  print ( '  Convection/diffusion equation with stabilization.' )
  print ( '  - ux - (mu+beta hmin) uxx = f in Omega = the unit interval.' )
  print ( '  u(0) = 0, u(1) = 1' )

  for my_grid in ( 10, 100 ):
    my_mu = 0.01
    my_beta = 0.5
    convection_diffusion_stabilized ( my_mu, my_grid, my_beta )
#
#  Terminate.
#
  print ( '' )
  print ( 'convection_diffusion_stabilized_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  convection_diffusion_stabilized_test ( )
