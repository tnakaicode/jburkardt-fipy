"""
Set up test #11.

  Discussion:

    Intersecting interfaces
    -del p del(u) = f on [-1,+1]x[-1,+1], where p is a piecewise constant function.
    u = g on the boundary.

    f = -del p del g
    g = r^alpha1 * mu(theta)

    Suggested parameter values:
      phi = pi/2
      p1 = R
      p2 = 1
      p3 = R
      p4 = 1 
      a1 = 0.1
      R = 161.4476387975881
      alpha1 = pi/4
      beta1 = -14.92256510455152.

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
theta = pi / 2
r = 161.4476387975881
alpha1 = pi / 4
beta1 = -14.92256510455152
p1 = r
p2 = 1.0
p3 = r
p4 = 1.0
a1 = 0.1
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
file = File ( "test11.pvd" )
file << u
#
#  Hold plot.
#
interactive ( )
