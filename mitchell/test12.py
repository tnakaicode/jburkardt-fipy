"""
Set up test #12.

  Multiple difficulties
  -Laplace(u) = f on the upside down L shaped region
  u = u0 on the boundary.
  u0 = awful formula
  f = -(u'')

  Discussion:

    Suggested parameter values:
      omega = 3 * pi / 2
      xw =  0.0
      yw = -0.75
      r0 =  0.75
      alphaw = 200.0
      xp = sqrt(5)/4
      yp = -0.25
      alphap = 1000
      epsilon = 0.01

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
omega = 3.0 * pi / 2.0
xw =  0.0
yw = -0.75
r0 =  0.75
alphaw = 200.0
xp = sqrt ( 5.0 ) / 4.0
yp = -0.25
alphap = 1000
epsilon = 0.01
#
#  Use a precomputed mesh.
#
mesh = Mesh ( "test12.xml" )
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
#
#  The right hand side function is simply -Laplacian(u0).
#
f = Expression ( "?" )
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
file = File ( "test12.pvd" )
file << u
#
#  Hold plot.
#
interactive ( )
