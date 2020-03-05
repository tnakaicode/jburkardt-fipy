"""
Set up test #4.

  Peak
  -Laplace(u) = f on the unit square.
  u = u0 on the boundary.
  u0 = exp(-alpha*(x-xc)^2+(y-yc)^2)
  f = -(u0'')

  Discussion:

    Suggested parameter values:
    * alpha = 1000,   (xc,yc) = (0.5,0.5)
    * alpha = 100000, (xc,yc) = (0.51,0.117)

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
#  Set parameters
#
alpha = 1000
xc = 0.5
yc = 0.5
#
#  Divide the unit square into a 10x10 mesh of quadrilaterals, and
#  split each into triangles.
#
mesh = UnitSquareMesh ( 10, 10 )
#
#  Define the function space.
#
V = FunctionSpace ( mesh, "Lagrange", 1 )
#
#  Define the exact solution.
#  Frap-dapping FENICS is failing to compile this expression.
#
u0 = Expression ( "exp ( - alpha * ( pow ( x[0] - xc, 2 ) + pow ( x[1] - yc, 2 ) ) )", \
                  alpha = alpha, \
                  xc = xc, \
                  yc = yc )

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
#  Define the right hand side, which is simply -(u0''):
#
f = Expression ( "4.0 * alpha * ( 1.0 - alpha "\
                 "* ( pow ( x[0] - xc, 2 ) + pow ( x[1] - yc, 2 ) ) )"\
                 "* exp ( - alpha * ( pow ( x[0] - xc, 2 ) + pow ( x[1] - yc, 2 ) ) )", \
                 alpha = alpha, \
                 xc = xc, \
                 yc = yc )

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
file = File ( "test04.pvd" )
file << u
#
#  Hold plot.
#
interactive ( )
