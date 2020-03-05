#! /usr/bin/env python
#
from fenics import *
from mshr import *

def fenics_to_txyz ( ):

#*****************************************************************************80
#
## fenics_to_txyz converts fenics data to nodes, elements, and values files.
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
# Define the variational problem.
#
  u = TrialFunction ( V )
  v = TestFunction ( V )
  f = Expression ( "10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree = 10 )
  g = Expression ( "sin(5*x[0])", degree = 10 )
  Auv = inner ( grad ( u ), grad ( v ) ) * dx
  Lv = f * v * dx + g * v * ds
#
# Compute the solution.
#
  u = Function ( V )
  solve ( Auv == Lv, u, bc )
#
#  Extract some information.
#
  u_nodal_values = u.vector()
  u_array = u_nodal_values.get_local()

  coords = mesh.coordinates()
  num_cells = mesh.num_cells()
  num_vertices = mesh.num_vertices()
  tri = mesh.cells()
#
#  Write the T file.
#
  tri_filename = 't.txt'

  outfile = open ( tri_filename,'w' )

  for row in tri:
    for column in row:
      outfile.write ( '%d\t' % column )
    outfile.write ( '\n' )

  outfile.close()

  print ( "Wrote triangle file '%s'" % ( tri_filename ) )
#
#  Write the XYZ file.  
#
  XYZ = np.zeros ( (num_vertices,3) )

  for k in range ( 0, num_vertices ):
    XYZ[k][0] = coords[k][0]
    XYZ[k][1] = coords[k][1]
    XYZ[k][2] = u_array[k]

  xyz_filename = 'xyz.txt'
  outfile = open ( xyz_filename, 'w' )

  for row in XYZ:
    for column in row:
      outfile.write ( '%14.8f' % column )
    outfile.write ( '\n' )

  outfile.close()

  print ( "Wrote xyz file '%s'" % ( xyz_filename ) )
#
#  Terminate.
#
  return

def fenics_to_txyz_test ( ):

#*****************************************************************************80
#
## fenics_to_txyz_test tests fenics_to_txyz.
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
  import time

  print ( time.ctime ( time.time() ) )
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  print ( '' )
  print ( 'fenics_to_txyz_test:' )
  print ( '  FENICS/Python version' )
  print ( "  Convert FENICS data to FEM files: triangles, nodes+values." )

  fenics_to_txyz ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'fenics_to_txyz_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  fenics_to_txyz_test ( )
