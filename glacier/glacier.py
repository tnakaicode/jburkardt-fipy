from dolfin import *
#
#  Discussion:
#
#    This FENICS input file sets up and solves the Stokes equations
#    for a model of a glacier.
#
#    The glacier is considered to lie within a rectangle tilted at
#    an angle ALPHA.
#
#  Modified:
#
#    11 February 2015
#
#  Reference:
#
#    William Mitchell,
#    Exact and numerical solutions for Stokes flows in glaciers,
#    Master's thesis,
#    Department of Mathematics and Statistics,
#    University of Alaska, Fairbanks, Alasaka,
#    August 2012
#
#  Set domain parameters and physical constants.
#
#  Length and Height.
#
Le = 4000.0
He = 500.0
#
#  Slope angle in radians.
#
alpha = pi / 180.0
#
#  Material density in kg/m^3.
#
rho = 917.0
#
#  Gravitational coefficient, m/s^2.
#
g = 9.81
#
#  Material viscosity, Pa-sec.
#
mu = 1.0E14
#
#  The body force vector G = (Gx,Gy).
#
G = Constant((sin(alpha) * g * rho, cos(alpha) * g * rho))
#
#  Define the mesh, using 3 cells in the X and Y direction.
#
mesh = RectangleMesh(0.0, 0.0, Le, He, 3, 3)
#
#  Specify the boundary conditions.
#


def LowerBoundary(x, on_boundary):
    return x[1] < DOLFIN_EPS and on_boundary


class PeriodicBoundary_x (SubDomain):
    def inside(self, x, on_boundary):
        return x[0] == 0 and on_boundary

    def map(self, x, y):
        y[0] = x[0] - Le
        y[1] = x[1]


pbc_x = PeriodicBoundary_x()

SlipRate = Expression(
    ("(3.0 + 1.7 * sin ( 2.0 * pi / %s * x[0] ) ) / 31557686.4" % Le, "0.0"))
#
#  Define the function spaces:
#    Velocity: piecewise quadratic vector.
#    Pressure: piecewise linear scalar.
#
V = VectorFunctionSpace(mesh, "CG", 2, constrained_domain=pbc_x)
Q = FunctionSpace(mesh, "CG", 1, constrained_domain=pbc_x)
W = V * Q
#
#  Define the Dirichlet condition at the base of the glacier.
#
bcD = DirichletBC(W.sub(0), SlipRate, LowerBoundary)
#
#  Define the periodic condition on the sides.
#

#bcP = PeriodicBC ( W.sub ( 0 ), pbc_x )
#
#  Define the variational problem: a(u,v) = L(v).
#
(v_i, q_i) = TestFunctions(W)
(u_i, p_i) = TrialFunctions(W)
a = (0.5 * mu * inner(grad(v_i) + grad(v_i).T, grad(u_i) + grad(u_i).T)
     - div(v_i) * pi + q_i * div(u_i)) * dx
L = inner(v_i, G) * dx
#
#  Matrix assembly.
#
U = Function(W)
#
#  Solution.
#
solve(a == L, U, bcD)
#
#  Split the mixed solution to recover velocity and pressure.
#
#( u, p ) = U.split ( )
# plot(u)
