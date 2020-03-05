"""
Set up test #6.

  Boundary layer
  -eps * Laplace(u) + 2 du/dx + du/dy = f on [-1,+1]x[-1,+1]
  u = u0 on the boundary.
  u0 = (1-exp(-(1-x)/eps)*(1-exp(-(1-y)/eps)*cos(pi(x+y))
  f = -eps*(u0'')+2 du0/dx+du0/dy

  Discussion:

    Suggested parameter values:
    * eps = 0.1
    * eps = 0.001

  Modified:

    30 August 2014

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
eps = 0.1
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
u0 = Expression ( "( 1.0 - exp ( - ( 1.0 - x[0] ) / eps ) ) * "
                  "( 1.0 - exp ( - ( 1.0 - x[1] ) / eps ) ) * "
                  "cos ( pi * ( x[0] + x[1] ) )" )
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
#  Define test and trial functions.
#
u = TrialFunction ( V )
v = TestFunction ( V )
#
#  Define the right hand side.
#
f = Expression 
  ( "exp ( - 2.0 / eps ) / pow ( eps, 2 ) * ("
    " - exp ( ( 1.0 - x[1] ) / eps )"
    " * ( -2.0 * exp ( 1.0 / eps ) * pow ( pi * eps, 2 )"
    " + exp ( x[0] / eps ) * ( 1.0 + 2.0 * pow ( pi * eps, 2 ) ) ) "
    " * cos ( pi * ( x[0] + x[1] ) ) "
    " - ( 3.0 * exp ( 2.0 / eps ) - exp ( ( 1.0 - x[0] ) / eps )"
    "  - exp ( ( 1.0 - x[1] ) / eps ) - exp ( ( x[0] + x[1] ) / eps ) )"
    "  * eps * pi * sin ( pi * ( x[0] + x[1] ) ) )",
    eps = eps,
    pi = pi
  )
#
#  Define the variational problem.
#  How can I BREAK UP THIS LONG LINE?
#  How do I express 'DUDX'?
#
a = ( eps * inner ( nabla_grad ( u ), nabla_grad ( v ) ) + 2.0 * inner ( DUDX, v ) + inner ( DUDY, v ) ) * dx
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
file = File ( "test06.pvd" )
file << u
#
#  Hold plot.
#
interactive ( )
