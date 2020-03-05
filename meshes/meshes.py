#! /usr/bin/env python
#
from fenics import *
from mshr import *

def meshes ( ):

#*****************************************************************************80
#
## meshes tries out various functions for meshing simple regions.
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
#  Unit interval [0,1]
#
  mesh = UnitIntervalMesh ( 10 )
  print ( 'Plotting a UnitIntervalMesh' )
  plot ( mesh, title = 'UnitIntervalMesh' )
  filename = 'unitintervalmesh.png'
  plt.savefig ( filename )
  print ( 'Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Any interval [a,b]
#
  a = -2.0
  b = 8.0
  mesh = IntervalMesh ( 10, a, b )
  print ( 'Plotting an IntervalMesh'  )
  plot ( mesh, title = 'IntervalMesh' )
  filename = 'intervalmesh.png'
  plt.savefig ( filename )
  print ( 'Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Unit square
#
  mesh = UnitSquareMesh ( 10, 10 )
  print ( 'Plotting a UnitSquareMesh' )
  plot ( mesh, title = 'UnitSquareMesh' )
  filename = 'unitsquaremesh.png'
  plt.savefig ( filename )
  print ( 'Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Unit square with "left" option.
#
  mesh = UnitSquareMesh ( 10, 10, 'left' )
  print ( 'Plotting a UnitSquareMesh' )
  plot ( mesh, title = 'UnitSquareMesh (left)' )
  filename = 'unitintervalmeshleft.png'
  plt.savefig ( filename )
  print ( 'Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Unit square with "crossed" option.
#
  mesh = UnitSquareMesh ( 10, 10, 'crossed' )
  print ( 'Plotting a UnitSquareMesh' )
  plot ( mesh, title = 'UnitSquareMesh (crossed)' )
  filename = 'unitintervalmeshcrossed.png'
  plt.savefig ( filename )
  print ( 'Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Unit square with "right/left" option.
#
  mesh = UnitSquareMesh ( 10, 10, 'right/left' )
  print ( 'Plotting a UnitSquareMesh' )
  plot ( mesh, title = 'UnitSquareMesh (right/left)' )
  filename = 'unitsquaremeshrightleft.png'
  plt.savefig ( filename )
  print ( 'Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Rectangle with SW and NE corners specified.
#
  sw = Point ( 0.0, 0.0 )
  ne = Point ( 10.0, 4.0 )
  mesh = RectangleMesh ( sw, ne, 10, 10 )
  print ( 'Plotting a RectangleMesh' )
  plot ( mesh, title = 'RectangleMesh')
  filename = 'rectanglemesh.png'
  plt.savefig ( filename )
  print ( 'Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Rectangle with "right/left" option.
#
  sw = Point ( -3.0, 2.0 )
  ne = Point ( 7.0, 6.0 )
  mesh = RectangleMesh ( sw, ne, 10, 10, 'right/left' )
  print ( 'Plotting a RectangleMesh' )
  plot ( mesh, title = 'RectangleMesh (right/left)' )
  filename = 'rectanglemeshrightleft.png'
  plt.savefig ( filename )
  print ( 'Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Circle.
#
  center = Point ( 1.0, 2.0 )
  mesh = generate_mesh ( Circle ( center, 1 ), 32 ) 
  print ( 'Plotting a Circle Mesh' )
  plot ( mesh, title = 'CircleMesh' )
  filename = 'circlemesh.png'
  plt.savefig ( filename )
  print ( 'Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Unit cube.
#
  mesh = UnitCubeMesh ( 10, 10, 10 )
  print ( 'Plotting a UnitCubeMesh' )
  plot ( mesh, title = 'UnitCubeMesh' )
  filename = 'unitcubemesh.png'
  plt.savefig ( filename )
  print ( 'Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Box, with SWF and NEB corner points.
#
  swf = Point ( 0.0, 0.0, 0.0 )
  neb = Point ( 10.0, 4.0, 2.0 )
  mesh = BoxMesh ( swf, neb, 10, 10, 10 )
  print ( 'Plotting a BoxMesh' )
  plot ( mesh, title = 'BoxMesh' )
  filename = 'boxmesh.png'
  plt.savefig ( filename )
  print ( 'Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Sphere.
#
  center = Point ( 1.0, 2.0, 3.0 )
  mesh = generate_mesh ( Sphere ( center, 1 ), 32 )
  print ( 'Plotting a Sphere Mesh' )
  plot ( mesh, title = 'Sphere Mesh' )
  filename = 'spheremesh.png'
  plt.savefig ( filename )
  print ( 'Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Terminate.
#
  return

def meshes_test ( ):

#*****************************************************************************80
#
## meshes_test tests meshes.
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
  print ( 'meshes_test:' )
  print ( '  FENICS/Python version' )
  print ( '  Demonstrate various built-in mesh types. )

  meshes ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'meshes_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  meshes_test ( )

