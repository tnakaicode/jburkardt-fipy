"""
Set up test #5.

  Battery
  -d/dx p dudx - d/dy q dudy = f on [0,8.4]x[0,24]
  Mixed boundary conditions
  u0 unknown
  f = given

  Discussion:

    P, Q and F are piecewise constant functions:

    Domain  P     Q       F
    ------ ----  ----    ----
    1      25.0  25.0     0.0
    2       7.0   0.8     1.0
    3       5.0   0.0001  1.0
    4       0.2   0.2     0.0
    5       0.05  0.05    0.0

   Y8  +------------------+
       | K = 1            |
   Y7  +---------------+  |
       | K = 5         |  |
   Y6  +--------+---+  |  |
       | K = 2  |   |  |  |
   Y5  +--------+ K |  |  |
       | K = 3  | = |  |  |
   Y4  +--------+ 4 |  |  |
       | K = 2  |   |  |  |
   Y3  +--------+   |  |  |
       | K = 5  |   |  |  |
   Y2  +--------+---+--+  |
       |                  |
   Y1  +------------------+
      X1       X2  X3 X4 X5

    The rectangular subdomains are characterized by 5 X values:
      0.0, 6.1, 6.5, 8.0, 8.4
    and 8 Y values:
      0.0, 0.8, 1.6, 3.6, 18.8, 21.2, 23.2, 24.0

    The boundary conditions are of the form:
      p dudx dyds - q dudy dxds + c * u = g

    using the following values:

      side     C    gN    Boundary condition
      ------- ---   ---   ------------------
      bottom  3.0   1.0:  - q dudy + 3 u = 1
      right   2.0   2.0:    p dudx + 2 u = 2
      top     1.0   3.0:    q dudy +   u = 3
      left    0.0   0.0:  - p dudx       = 0

  Modified:

    23 August 2014

  Author:

    John Burkardt

  Reference:

    William Mitchell,
    A collection of 2D elliptic problems for testing adaptive
    grid refinement algorithms,
    Applied Mathematics and Computation,
    Volume 220, 1 September 2013, pages 350-364.
"""

from dolfin import *
#
#  Get the mesh from a file.
#  I had to create the mesh in MESHFACES,
#  write the nodes and elements to files,
#  then use TRIANGULATION_TO_XML to convert the nodes and elements to XML.
#
mesh = Mesh ( "mesh05.xml" )
#
#  Define the function space.
#
V = FunctionSpace ( mesh, "Lagrange", 1 )
#
#  Apply the boundary condition to points which are on the boundary of the mesh.
#
def u0_boundary ( x, on_boundary ):
  return on_boundary
#
#  REVISE THE BC HERE
#
bc = DirichletBC ( V, u0, u0_boundary )
#
#  Define the variational problem.
#
u = TrialFunction ( V )
v = TestFunction ( V )
#
#  THE RIGHT HAND SIDE IS PIECEWISE...
#
f = Expression ( "?" )
#
#  THE OPERATOR IS DIFFERENT.
#
a = inner ( nabla_grad ( u ), nabla_grad ( v ) ) * dx
L = f * v * dx
#
#  Compute the solution.
#
u = Function ( V )
solve ( a == L, u, bc )
#
#  Plot the solution and the mesh.
#
plot ( u )
plot ( mesh )
#
#  Dump solution to file in VTK format.
#  It can be viewed by programs such as ParaView or VisIt.
#
file = File ( "test05.pvd" )
file << u
#
#  Hold plot.
#
interactive ( )
