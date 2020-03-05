#! /usr/bin/env python
#
from fenics import *

def expression_test1 ( ):

#*****************************************************************************80
#
## expression_test1 demonstrates some features of the Expression() function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    06 November 2018
#
#  Author:
#
#    John Burkardt
#
  import matplotlib.pyplot as plt

  print ( '' )
  print ( 'expression_test1:' )
#
#  Define an expression using "degree".
#
  f_expr = Expression ( 'sin ( pi * x[0] ) * sin ( 2 * pi * x[1] )', degree = 10 )
#
#  We can ask for the type of f_expr.
#
  print ( '' )
  print ( '  Request the type of the expression object:' )
  print ( '  type ( f_expr ) = ', type ( f_expr ) )
#
#  The "print" command does not return the definition string.
#  I don't know how to get it back.
#
  print ( '' )
  print ( '  The print command does not return the definition string.' )
  print ( '  print ( f_expr ) = ', f_expr )
#
#  Evaluate the expression.
#
  print ( '' )
  print ( '  Choose a single argument for the function:' )
  x = ( 0.40, 0.24 )
  print ( '  x = ', x )

  print ( '' )
  print ( '  Evaluate the function:' )
  value = f_expr ( x )
  print ( '  f_expr ( x ) = ', value )
#
#  Plot the expression over a given mesh.
#
  my_mesh = UnitSquareMesh ( 10, 10 )
  fig = plot ( f_expr, mesh = my_mesh, title = 'expression_test1' )
  filename = 'expression_test1.png'
  plt.savefig ( filename )
  print ( '  Saving graphics in "%s"' % filename )
#
#  Create an expression with default parameter values:
#
  print ( '' )
  print ( '  Define g_expr = Expression ( "pow ( x[0], POWER )", POWER = 2, ... )' )
  g_expr = Expression ( "pow ( x[0], POWER )", POWER = 2, degree = 10 )
#
#  Evaluate the expression with the default parameter value.
#
  x = 3.0
  print ( '  x = ', x )
  value = g_expr ( x )
  print ( '  g_expr ( x ) = ', value )
#
#  Change the parameter value and reevaluate.
#
  print ( '  Now reset g_expr.POWER = 3' )
  g_expr.POWER = 3
  value = g_expr ( x )
  print ( '  g_expr ( x ) = ', value )
#
#  Define a vector-valued expression.
#
  print ( '' )
  print ( '  Define h_expr = Expression ( ( "sin(pi*x[0])", "cos(pi*x[0])" ), degree = 10 )' )

  h_expr = Expression ( ( 'sin(pi*x[0])', 'cos(pi*x[0])' ), degree = 10 )
  x = 1.0 / 3.0
  value = h_expr ( x )
  print ( '  h_expr ( x ) = ', value )
#
#  Terminate.
#
  return

def expression_test2 ( ):

#*****************************************************************************80
#
## expression_test2 demonstrates the use of Expression() function.
#
#  Discussion:
#
#    We set up the following Poisson problem:
#
#    - div grad u = - 6         in Omega
#    u(dOmega) = u_exact(x,y)   on dOmega
#
#    u_exact = 1 + x^2 + 2 y^2
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
#    Hans Petter Langtangen, Anders Logg,
#    Solving PDEs in Python - The FEniCS Tutorial Volume !
#
  import matplotlib.pyplot as plt

  print ( '' )
  print ( 'expression_test2:' )
  print ( '  Use Expression() in a finite element calculation.' )
  print ( '' )
#
#  Create an 8x8 triangular mesh on the unit square.
#
  my_mesh = UnitSquareMesh ( 8, 8 )
#
#  Define the function space.
#
  V = FunctionSpace ( my_mesh, 'P', 1 )
#
#  Define the exact solution using an expression:
#
  u_exact_expr = Expression ( '1 + x[0]*x[0] + 6*x[1]*x[1]', degree = 4 )
#
#  Use the expression to define the boundary condition.
#
  def boundary ( x, on_boundary ):
    return on_boundary

  bc = DirichletBC ( V, u_exact_expr, boundary )
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
  uh = Function ( V )
  solve ( a == L, uh, bc )
#
#  Use the expression to compute the error in the L2 norm.
#
  error_L2 = errornorm ( u_exact_expr, uh, 'L2' )
  print ( '  error_L2  =', error_L2 )
#
#  Use the expression to compute the maximum error at mesh vertices.
#
  vertex_values_u_exact_expr = u_exact_expr.compute_vertex_values ( my_mesh )
  vertex_values_uh = uh.compute_vertex_values ( my_mesh )
  import numpy as np
  error_max = np.max ( np.abs ( vertex_values_u_exact_expr - vertex_values_uh ) )
  print ( '  error_max =', error_max )
#
#  Terminate.
#
  return

def expression_test3 ( ):

#*****************************************************************************80
#
## expression_test3 demonstrates subclassing an Expression.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    06 November 2018
#
#  Author:
#
#    John Burkardt
#
  import matplotlib.pyplot as plt

  print ( '' )
  print ( 'expression_test3:' )
  print ( '  Use subclassing to set up an expression that is too complicated' )
  print ( '  for a one-line definition.' )
#
#  Mesh the unit square.
#
  my_grid = 4
  my_mesh = UnitSquareMesh ( my_grid, my_grid, 'right/left' )
#
#  Define the function space.
#
  V = FunctionSpace ( my_mesh, 'P', 1 )
#
#  Note that we can use the Python exponentiation format here.
#  I'm not sure where this is actually explained, with a comprehensible example.
#
  class k_user_expr ( UserExpression ):
    def eval ( self, values, x ):
      t1 = - 8 * ( pi * ( x[0]**2 - 2 * x[0] ) )**2 + cos ( 4 *pi*x[1] )
      t2 = ( x[0]**4 - 4*x[0]**3 - 8*x[0]**2 + 24*x[0] - 8 ) * sin(2*pi*x[1])**2
      values[0] = t1 + t2
    def value_shape ( self ):
      return ()
#
#  Turn the user expression into an expression by including the degree.
#
  k_expr = k_user_expr ( degree = 4 )
#
#  Now we can evaluate the expression, if we include the degree.
#
  x = ( 0.40, 0.24 )
  print ( '  x = ', x )
  k_value = k_user_expr ( x, degree = 4 )
  print ( '  k_expr ( x ) = ', k_value )
#
#  Project the expression into V, making a finite element function.
#
  k_fefun = project ( k_user_expr ( ), V )
#
#  Plot the finite element function.
#
  fig = plot ( k_fefun, title = 'k_fefun' )
  plt.colorbar ( fig )
  filename = 'expression_test3.png'
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Terminate.
#
  return

def expression_test ( ):

#*****************************************************************************80
#
## expression_test tests the FEniCS "Expression" function.
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
  import time
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  print ( time.ctime ( time.time() ) )
  print ( '' )
  print ( 'expression_test:' )
  print ( '  FENICS/Python version' )
  print ( '  Demonstrate the FEniCS "Expression()" function.' )

  expression_test1 ( )
  expression_test2 ( )
  expression_test3 ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'expression_test:' );
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )

if ( __name__ == '__main__' ):

  expression_test ( )
