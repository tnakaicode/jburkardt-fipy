"""
Set up test #2.

  Reentrant corner
  -Laplace(u) = f on [-1,+1]x[-1,+1] with a section removed from the clockwise
  side of the positive X axis with angle omega,
  u = u0 on the boundary.
  u0 = r^alpha sin(alpha*theta) 
  f = -(u0'')

  Set omega = 3pi/2 for the L-domain problem.

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
#  Define parameters.
#  The symbol "pi" is predefined.
#
omega = 3 * pi / 2
alpha = pi / omega
#
#  The precomputed mesh should be read from the given file.
#
mesh = Mesh ( "test02.xml" )
#
#  Define the function space.
#
V = FunctionSpace ( mesh, "Lagrange", 1 )
#
#  Define the exact solution.
#  I assume the function "atan2()" is available.
#
u0 = Expression ( "pow ( sqrt(x[0]*x[0]+x[1]*x[1]), alpha )" 
                  "* sin ( alpha * atan2(x[1],x[0]) )" )
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
#
#  Got to figure out the correct value of F here!
#
f = Expression ( "???" )

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
file = File ( "test02.pvd" )
file << u
#
#  Hold plot.
#
interactive ( )
