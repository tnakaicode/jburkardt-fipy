"""
Set up test #1.

  Analytic solution
  -Laplace(u) = f on the unit square.
  u = u0 on the boundary.
  u0 = 2^(4p)*x^p*(1-x)^p*y^p*(1-y)^p, 
  f = 0.

  Discussion:

    The parameter P should be chosen large enough so that the highest order
    finite elements to be used will not give the exact solution.  A value
    of P = 10 might be suitable.

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
#  Set parameter values.
#
p = 10
#
#  Divide the unit square [0,1]x[0,1] into a mesh of quadrilaterals, and
#  split each into triangles.
#
mesh = UnitSquareMesh ( 10, 10 )
#
#  Define the function space.
#
V = FunctionSpace ( mesh, "Lagrange", 1 )
#
#  Define the exact solution.
#
#  A long expression can be written on multiple lines, with quotes
#  surrounding each part.
$
#  If you use parameters in the expression, (such as P), then you
#  need to supply the value (as in "p = p"), but you must use a comma
#  to terminate the expression, and no quotes around the parameter line.
#
u0 = Expression ( 'pow ( 2, 4*p) * pow ( x[0], p ) * pow ( 1-x[0], p )'
                  '              * pow ( x[1], p ) * pow ( 1-x[1], p )',
                  p = p )
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
f = Expression ( "-(                                                "\
                 "       (  p*(p-1)*pow(x[0],p-2)*pow(1-x[0],p)     "\
                 "         -p* p   *pow(x[0],p-1)*pow(1-x[0],p-1)   "\
                 "         +p*(p-1)*pow(x[0],p)  *pow(1-x[0],p-2) ) "\
                 "     * pow(x[1],p)*pow(1-x[1],p)                  "\
                 "   +   pow(x[0],p)*pow(1-x[0],p)                  "\
                 "     * (  p*(p-1)*pow(x[1],p-2)*pow(1-x[1],p)     "\
                 "        - p* p   *pow(x[1],p-1)*pow(1-x[1],p-1)   "\
                 "        + p*(p-1)*pow(x[1],p)  *pow(1-x[1],p-2) ) )",\
   p = p )
#
a = inner ( nabla_grad ( u ), nabla_grad ( v ) ) * dx
L = f * v * dx
#
#  Specify that the solution u is to be a member of the function space V.
#
u = Function ( V )
#
#  Compute the solution u which solves the equation with given boundary conditions.
#
solve ( a == L, u, bc )
#
#  Plot the solution to the screen.
#
plot ( u )
#
#  Plot the mesh to the screen.
#
plot ( mesh )
#
#  Dump the solution to file in VTK format.
#  It can be viewed by programs such as ParaView or VisIt.
#
file = File ( "test01.pvd" )
file << u
#
#  Hold plot.
#
interactive ( )
