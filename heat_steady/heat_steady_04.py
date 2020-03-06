#! /usr/bin/env python3
#
from fenics import *


def heat_steady_04():

    #*****************************************************************************80
    #
    ## heat_steady_04, 2D steady heat equation on a rectangle, discontinuous diffusivity k.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    23 October 2018
    #
    #  Author:
    #
    #    John Burkardt
    #
    import matplotlib.pyplot as plt
    #
    #  Define the mesh.
    #
    x_left = 0.0
    x_right = 5.0
    y_bottom = 0.0
    y_top = 1.0

    sw = Point(x_left, y_bottom)
    ne = Point(x_right, y_top)

    mesh = RectangleMesh(sw, ne, 50, 10)
    #
    #  Define the function space.
    #
    V = FunctionSpace(mesh, "Lagrange", 1)
    #
    #  Define boundary conditions.
    #
    u_left = 100.0

    def on_left(x, on_boundary):
        return (x[0] <= x_left + DOLFIN_EPS)

    bc_left = DirichletBC(V, u_left, on_left)

    u_right = 0.0

    def on_right(x, on_boundary):
        return (on_boundary and x_right - DOLFIN_EPS <= x[0])

    bc_right = DirichletBC(V, u_right, on_right)

    bc = [bc_left, bc_right]
    #
    #  Define the trial functions (u) and test functions (v).
    #
    u = TrialFunction(V)
    v = TestFunction(V)
    #
    #  The diffusivity is a piecewise function defined by this
    #  ugly C expression.
    #
    k = Expression('x[0] < 2.0 ? 5.0 : 1.0', degree=0)
    #
    #  Write the bilinear form.
    #
    Auv = k * inner(grad(u), grad(v)) * dx
    #
    #  Write the linear form.
    #
    f = Constant(0.0)
    Lv = f * v * dx
    #
    #  Solve the variational problem with boundary conditions.
    #
    u = Function(V)
    solve(Auv == Lv, u, bc)
    #
    #  Plot the solution.
    #
    plot(u, title='heat_steady_04')
    filename = 'heat_steady_04.png'
    plt.savefig(filename)
    print("Saving graphics in file '%s'" % (filename))
    plt.close()
    #
    #  Terminate.
    #
    return


def heat_steady_04_test():

    #*****************************************************************************80
    #
    ## heat_steady_04_test tests ell.
    #
    #  Modified:
    #
    #    21 October 2018
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
    print('heat_steady_04_test:')
    print('  FENICS/Python version')
    print('  2D steady heat equation on a rectangle.')
    print('  Discontinuous diffusivity coefficient K.')

    heat_steady_04()
    #
    #  Terminate.
    #
    print('')
    print('heat_steady_04_test:')
    print('  Normal end of execution.')
    print('')
    print(time.ctime(time.time()))
    return


if (__name__ == '__main__'):

    heat_steady_04_test()
