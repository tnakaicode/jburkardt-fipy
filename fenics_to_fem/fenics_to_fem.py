#! /usr/bin/env python
#
from fenics import *
from mshr import *

def fenics_to_fem ( ):

#*****************************************************************************80
#
## fenics_to_fem converts fenics data to nodes, elements, and values files.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    17 October 2018
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import math 
#
#  Create a mesh.
#
  radius = 1.0
  domain = Circle ( Point ( 0.0, 0.0 ), radius, 12 )
  mesh = generate_mesh ( domain, 20 )
#
#  Define the function space.
#
  V = FunctionSpace ( mesh, "Lagrange", 1 )
#
#  Define the Dirichlet boundary.
#
  def boundary ( x, on_boundary ):
    return on_boundary
#
#  Define the boundary condition.
#
  u0 = Constant ( 0.0 )
  bc = DirichletBC ( V, u0, boundary )
#
#  Define the variational problem.
#
  u = TrialFunction ( V )
  v = TestFunction ( V )
  f = Expression ( "10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree = 10 )
  g = Expression ( "sin(5*x[0])", degree = 10 )
  Auv = inner ( grad ( u ), grad ( v ) ) * dx
  Lv = f * v * dx + g * v * ds
#
#  Compute the solution.
#
  u = Function ( V )
  solve ( Auv == Lv, u, bc )
#
#  Write NODES, ELEMENTS and VALUES to files.
#
  header = 'spike'
#
#  Create a NODES array, and write it to a file.
#
  node_num = mesh.num_vertices ( )
  node_xy = mesh.coordinates ( )

  node_filename = header + '_nodes.txt'
  node_file = open ( node_filename, 'w' )

  for row in node_xy:
    for column in row:
      node_file.write ( '  %g' % column )
    node_file.write ( '\n' )

  node_file.close ( )

  print ( "Wrote node file '%s'" % ( node_filename ) )
#
#  Create an ELEMENTS array, and write it to a file.
#
  element_num = mesh.num_cells ( )
  element_node = mesh.cells ( )
  element_filename = header + '_elements.txt'
  element_file = open ( element_filename, 'w' )

  for row in element_node:
    for column in row:
      element_file.write ( '  %d' % column )
    element_file.write ( '\n' )

  element_file.close ( )

  print ( "Wrote element file '%s'" % ( element_filename ) )
#
#  Create a VALUES array, and write it to a file.
#
  u_nodal_values = u.vector ( )
  u_array = u_nodal_values.get_local ( )

  value_num = mesh.num_vertices ( )
  node_value = np.zeros ( node_num )
  value_filename = header + '_values.txt'
  value_file = open ( value_filename, 'w' )

  for i in range ( 0, value_num ):
    value_file.write ( '  %g\n' % u_array[i] )

  value_file.close ( )

  print ( "Wrote value file '%s'" % ( value_filename ) )
#
#  Terminate.
#
  return

def fenics_to_fem_test ( ):

#*****************************************************************************80
#
## fenics_to_fem_test tests fenics_to_fem.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
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
  print ( 'fenics_to_fem_test:' )
  print ( '  FENICS/Python version' )
  print ( "  Convert FENICS data to FEM files: nodes, elements, values." )

  fenics_to_fem ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'fenics_to_fem_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  fenics_to_fem_test ( )
