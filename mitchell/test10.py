"""
Set up test #10.

  Interior line singularity
  -Laplace(u) = f on [-1,+1]x[-1,+1]
  u = u0 on the boundary.
  u0 = cos(pi*y/2)                        for x <= beta*(y+1)
     = cos(pi*y/2)+(x-beta*(y+1))**alpha  for      beta*(y+1) < x
  f = -Laplace(u).

  Discussion:

    Suggested parameter values:
    * alpha = 2.5, beta = 0.0
    * alpha = 1.1, beta = 0.0
    * alpha = 1.5, beta = 0.6

  Modified:

    22 August 2014

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
#  Set parameters.
#
alpha = 2.5
beta = 0.0
#
#  Define a mesh of triangles on [-1,+1]x[-1,+1]
#
mesh = RectangleMesh ( -1.0, -1.0, +1.0, +1.0, 4, 4 )
#
#  Define the function space.
#
V = FunctionSpace ( mesh, "Lagrange", 1 )
#
#  Define the exact solution.
#
u0 = Expression ( "?" )
#
#  Apply the boundary condition to points which are on the boundary of the mesh.
#
def u0_boundary ( x, on_boundary ):
  return on_boundary
#
#  The value to be used at Dirichlet boundary points is U0.
#
bc = DirichletBC ( V, u0, u0_boundary )
#
#  Define the variational problem.
#
u = TrialFunction ( V )
v = TestFunction ( V )
f = Expression ( "?" )
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
file = File ( "test10.pvd" )
file << u
#
#  Hold plot.
#
interactive ( )
