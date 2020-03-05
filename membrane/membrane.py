#! /usr/bin/env python
#
from fenics import *
from mshr import *


def membrane():

    # *****************************************************************************80
    #
    # membrane solves the Poisson equation on a circular membrane.
    #
    #  Modified:
    #
    #    23 October 2018
    #
    #  Author:
    #
    #    John Burkardt
    #
    #  Reference:
    #
    #    Hans Petter Langtangen, Anders Logg,
    #    Solving PDEs in Python - The FEniCS Tutorial Volume !
    #
    import matplotlib.pyplot as plt
    import numpy
#
#  Set the pressure function:
#  T = Tension,
#  A = Amplitude,
#  R = Radius of domain.
#
    T = 10.0
    A = 1.0
    R = 0.3
    theta = 0.2
    x0 = 0.6 * R * cos(theta)
    y0 = 0.6 * R * sin(theta)
#
#  Set sigma to a large value like 50 for verification.
#
    sigma = 0.025
#
#  Set the approximate number of elements in the radial direction.
#
    n = 40
    center = Point(0.0, 0.0)
    mesh = generate_mesh(Circle(center, 1.0), n)
#
#  Define the function space to be used for Trial, Test, and Solution.
#
    V = FunctionSpace(mesh, "Lagrange", 1)
#
#  Define boundary condition w = 0.
#

    def boundary(x, on_boundary):
        return on_boundary

    bc = DirichletBC(V, Constant(0.0), boundary)
#
#  Define trial and test functions.
#
    w = TrialFunction(V)
    v = TestFunction(V)
#
#  Define the variational problem.
#
    auv = inner(nabla_grad(w), nabla_grad(v)) * dx

    f = Expression(
        "4*exp(-0.5*(pow((R*x[0] - x0)/sigma, 2 )) "
        "      -0.5*(pow((R*x[1] - y0)/sigma, 2 ))) ",
        R=R, x0=x0, y0=y0, sigma=sigma, degree=10)

    fv = f * v * dx
#
#  Compute the solution.
#
    w = Function(V)

    problem = LinearVariationalProblem(auv, fv, w, bc)

    solver = LinearVariationalSolver(problem)
    solver.parameters["linear_solver"] = "cg"
    solver.parameters["preconditioner"] = "ilu"
    solver.solve()
#
#  Plot the mesh.
#
    plot(mesh, title="Mesh over domain")
    filename = 'membrane_mesh.png'
    plt.savefig(filename)
    print('Graphics saved as "%s"' % (filename))
    plt.close()
#
#  Plot the solution.
#
    plot(w, title="Scaled deflection")
    filename = 'membrane_solution.png'
    plt.savefig(filename)
    print('Graphics saved as "%s"' % (filename))
    plt.close()
#
#  Plot the pressure.
#
    f = interpolate(f, V)
    plot(f, title="Scaled pressure")
    filename = 'membrane_pressure.png'
    plt.savefig(filename)
    print('Graphics saved as "%s"' % (filename))
    plt.close()
#
#  Find the maximum real deflection.
#
    max_w = w.vector().get_local().max()
    max_D = A * max_w / (8 * pi * sigma * T)
    print('')
    print("  Maximum real deflection is %g" % (max_D))
#
#  Verification for "flat" pressure (large sigma).
#
    if 50 <= sigma:
        w_exact = Expression("1 - x[0]*x[0] - x[1]*x[1]")
        w_e = interpolate(w_exact, V)
        dev = numpy.abs(w_e, vector().get_local() -
                        w.vector().get_local()).max()
        print('sigma=%g: max deviation=%e' % (dev))
#
#  Terminate.
#
    return


def membrane_test():

    # *****************************************************************************80
    #
    # membrane_test tests membrane.
    #
    #  Modified:
    #
    #    23 October 2018
    #
    #  Author:
    #
    #    John Burkardt
    #
    import time

    print(time.ctime(time.time()))
#
#  Report level = only warnings or higher.
#
    level = 30
    set_log_level(level)

    print('')
    print('membrane_test:')
    print('  FENICS/Python version')
    print('  Poisson equation on a circular membrane.')

    membrane()
#
#  Terminate.
#
    print('')
    print('membrane_test:')
    print('  Normal end of execution.')
    print('')
    print(time.ctime(time.time()))
    return


if (__name__ == '__main__'):

    membrane_test()
