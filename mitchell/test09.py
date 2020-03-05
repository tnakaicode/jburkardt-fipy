"""
Set up test #9.

  Wave front
  -Laplace(u) = f on the unit square.
  u = u0 on the boundary.
  u0 = atan ( alpha * (sqrt(x-xc)**2+(y-yc)**2)-r0) )
  f = -Laplace(u0)

  Discussion:

    Suggested parameter values:
    * alpha =   20, xc = -0.05, yc = -0.05, r0 = 0.7
    * alpha = 1000, xc = -0.05, yc = -0.05, r0 = 0.7
    * alpha = 1000, xc =  1.5,  yc =  0.25, r0 = 0.92
    * alpha =   50, xc =  0.5,  yc =  0.5,  r0 = 0.25

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
#  Set parameters.
#
alpha = 20.0
xc = -0.05
yc = -0.05
r0 = 0.7
#
#  Divide the unit square into a 4x4 mesh of quadrilaterals, and
#  split each into triangles.
#
mesh = UnitSquareMesh ( 4, 4 )
#
#  Define the function space.
#
V = FunctionSpace ( mesh, "Lagrange", 1 )
#
#  Define the radial coordinate R.
#
r = Expression ( "sqrt ( pow ( x[0] - xc, 2 ) + pow ( x[1] - yc, 2 ) - r0 )" ,
                 xc = xc, \
                 yc = yc, \
                 r0 = r0 )
#
#  Define the exact solution.
#  Can we use r in this way?  Do we need "r = r"?
#
u0 = Expression ( " atan ( alpha * r )" ,
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
#  Define F.
#
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
file = File ( "test09.pvd" )
file << u
#
#  Hold plot.
#
interactive ( )
