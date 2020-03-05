"""
Set up test #8.

  Oscillatory
  -Laplace(u) -u/(alpha+r)^4 = f on the unit square.
  r = sqrt ( x^2 + y^2 )
  u = u0 on the boundary.
  u0 = sin(1/(alpha+r))
  f = -Laplace(u0) -u0/(alpha+r)^4

  Discussion:

    The parameter alpha affects the number of oscillations in the solution.

    Suggested parameter values:
    * alpha = 1/(10*pi) is relatively easy.
    * alpha = 1/(50*pi) is harder.

  Modified:

    27 August 2014

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
alpha = 1.0 / ( 10.0 * pi )
#
#  Divide the unit square [0,1]x[0,1] into a mesh of quadrilaterals, and
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
u0 = Expression ( "sin ( 1.0 / ( alpha + sqrt ( pow ( x[0], 2 ) + pow ( x[1], 2 ) ) )",
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
#  Define the radial coordinate R.
#
r = Expression ( "sqrt ( x[0] * x[0] + x[1] * x[1] )" )
#
#  The right hand side:
#
f = Expression ( "?, "\                                               "\
   "alpha = alpha" )
#
#  Define the operator.
#  Do we use **4 here?  What blasted language are we talking?
#
a = ( inner ( nabla_grad ( u ), nabla_grad ( v ) ) - u / ( alpha + r )**4 ) * dx
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
file = File ( "test08.pvd" )
file << u
#
#  Hold plot.
#
interactive ( )
