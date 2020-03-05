"""
Set up test #3.

  Linear Elasticity
  Coupled system of equations, to be described later.
  Region is the square [-1,+1]x[-1,+1] with a slit from (0,0) to (1,0).
  u = u0 on the boundary
  u0 ...
  f ...

  Suggested Parameter Values:

    lambda = 0.5444837367825, q =  0.5430755788367
    lambda = 0.9085291898461, q = -0.2189232362488

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
lambda = 0.5444837367825
q = 0.5430755788367
#
#  Read the precomputed mesh.
#
mesh = Mesh ( "test03.xml" )
#
#  Define the function space.
#
V = FunctionSpace ( mesh, "Lagrange", 1 )
#
#  Define the exact solution.
#  HOW DO YOU DEAL WITH U0 BEING A VECTOR?
#  OR DO WE SET U0 AND V0?
#
u0 = Expression ( "?" )
#
#  NEED TWO COMPONENTS HERE
#
def u0_boundary ( x, on_boundary ):
  return on_boundary
#
#  ALL THE FOLLOWING MUST BE MODIFIED FOR TWO COMPONENT PROBLEM
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
file = File ( "test03.pvd" )
file << u
#
#  Hold plot.
#
interactive ( )
