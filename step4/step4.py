#! /usr/bin/env python
#
from fenics import *

def step4 ( ):

#*****************************************************************************80
#
## step4 shows how to generate the right hand side of a Poisson problem.
#
#  Discussion:
#
#    We are interested in solving an equation of the form:
#
#      - div k(u,x,y) grad u = f in Omega
#                          u = g on dOmega
#
#    We have chosen:
#      u(x,y) = 1 + x + 2 * y
#      k(u,x,y) = 1 + u^2
#      
#    The corresponding f(x,y) is now determined.  But instead of computing it
#    by hand, we ask the Python symbolic package to discover the formula:
#      f(x,y) = - 10 * x - 20 * y - 10
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    25 October 2018
#
#  Author:
#
#    John Burkardt
#
  import matplotlib.pyplot as plt
  import sympy as sym

  def k ( u ):
    value = 1 + u ** 2
    return value

  x, y = sym.symbols ( 'x[0], x[1]' )
  u = 1 + x + 2 * y

  f = - sym.diff ( k ( u ) * sym.diff ( u, x ), x ) \
      - sym.diff ( k ( u ) * sym.diff ( u, y ), y )

  f = sym.simplify ( f )

  f_code = sym.printing.ccode ( f )
  print ( '' )
  print ( '  Formula for F(x,y) computed by SymPy:' )
  print ( '' )
  print ( '  f(x,y) = ', f_code )

  u_code = sym.printing.ccode ( u )

#----------------------------------
#  Can now set up and solve problem!
#-----------------------------------
#
#  Create an 8x8 triangular mesh on the unit square.
#
  mesh = UnitSquareMesh ( 8, 8 )
#
#  Define the function space.
#
  V = FunctionSpace ( mesh, 'P', 1 )
#
#  Define the exact solution using the form generated above.
#
  u_D = Expression ( u_code, degree = 1 )
#
#  Define the boundary condition, using the exact solution.
#
  def boundary ( x, on_boundary ):
    return on_boundary

  bc = DirichletBC ( V, u_D, boundary )
#
#  Define the variational problem, using the formula for f.
#
  u = Function ( V )
  v = TestFunction ( V )
  f = Expression ( f_code, degree = 1 )
  F = k ( u ) * dot ( grad ( u ), grad ( v ) ) * dx - f * v * dx
#
#  Compute the solution.
#
#  Because we no longer have a linear system (that is, k(u) depends on u), 
#  we have to ask FENICS to solve F==0, not a==L!
#
  solve ( F == 0, u, bc )
#
#  Plot the solution.
#
  plot ( u, mode = 'contour', title = 'step4 solution' )
  filename = 'step4_solution.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Terminate.
#
  return

def step4_test ( ):

#*****************************************************************************80
#
## step4_test tests step4.
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
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  print ( '' )
  print ( 'step4_test:' )
  print ( '  FENICS/Python version' )
  print ( '  Generate the right hand side function f(x,y)' )
  print ( '  for a Poisson problem' )
  print ( '  - div k(u,x,y) grad u = f(x,y) in Omega' )
  print ( '                      u = g(x,y) on dOmega' )
  print ( '  where we have already chosen' )
  print ( '    u(x,y) = 1 + x + 2y' )
  print ( '    k(u,x,y) = 1 + u^2' )
  print ( '  and want f(x,y) computed symbolically by SymPy.' )

  step4 ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'step4_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  step4_test ( )
