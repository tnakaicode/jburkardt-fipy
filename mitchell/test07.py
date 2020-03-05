"""
Set up test #7.

  Boundary line singularity
  -Laplace(u) = f on the unit square.
  u = u0 on the boundary.
  u0 = x^alpha
  f = - alpha * ( alpha - 1 ) * x^(alpha-2)

  Discussion:

    The parameter alpha should be 0.5 or greater.  It determines the
    strength of the singularity.

    Suggested parameter values:
    * alpha = 0.6

  Modified:

    24 August 2014

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
alpha = 0.6
#
#  Divide the unit square [0,1]x[0,1] into a 4x4 mesh of quadrilaterals, and
#  split each into triangles.
#
mesh = UnitSquareMesh ( 4, 4 )
#
#  Define the function space.
#
V = FunctionSpace ( mesh, "Lagrange", 1 )
#
#  Define the exact solution.
#
u0 = Expression ( "pow ( x[0], alpha )",
                  alpha = alpha )
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
#  We have a singularity in F at X = 0.
#
f = Expression ( "- alpha * ( alpha - 1.0 ) * pow ( x[0]+0.000001, alpha - 2 )",
                 alpha = alpha )
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
file = File ( "test07.pvd" )
file << u
#
#  Hold plot.
#
interactive ( )
