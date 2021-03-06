#! /usr/bin/env python
#
from fenics import *


def bvp_04(e_num=8):

    # *****************************************************************************80
    #
    # bvp_04 solves a boundary value problem.
    #
    #  Discussion:
    #
    #    -u'' + u = x, 0 < x < 1
    #    u(0) = 0, u(1) = 0
    #
    #    Exact solution is u(x) = x - sinh(x)/sinh(1)
    #
    #    This function solves the same problem as bvp_01.py; however,
    #    the Dirichlet boundary conditions are now specified separately
    #    and explicitly.  This is done in order that we can next show how
    #    to specify a Neumann condition at one endpoint.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    17 January 2016
    #
    #  Author:
    #
    #    Hans-Werner van Wyk
    #
    import matplotlib.pyplot as plt
#
#  Define the mesh
#
    mesh = UnitIntervalMesh(e_num)
#
#  Function Space
#
    V = FunctionSpace(mesh, "CG", 1)
#
#  Define the right hand side.
#  (The "domain = mesh" argument stops FENICS from complaining that it
#  isn't sure what the spatial dimension is...)
#
    f = Expression("x[0]", degree=10)
#
#  Define the exact Function.
#
    ue = Expression("x[0] - sinh ( x[0] ) / sinh ( 1.0 )", degree=1)
#
#  Left boundary.
#

    def left_boundary(x, on_boundary):
        return on_boundary and (abs(x[0] - 0.0) < DOLFIN_EPS)
#
#  Right boundary.
#

    def right_boundary(x, on_boundary):
        return on_boundary and (abs(x[0] - 1.0) < DOLFIN_EPS)

    x_left = 0.0
    u_left = ue(x_left)
    x_right = 1.0
    u_right = ue(x_right)

    bc1 = DirichletBC(V, u_left, left_boundary)
    bc2 = DirichletBC(V, u_right, right_boundary)

    bc = [bc1, bc2]
#
#  Define the weak form
#
    psii = TrialFunction(V)
    psij = TestFunction(V)

    A = inner(grad(psii), grad(psij)) * dx + psii * psij * dx
    RHS = f * psij * dx
#
#  Solve
#
    u = Function(V)

    solve(A == RHS, u, bc)
#
#  Plot the solution.
#
    plot(u, title='bvp_04 solution')
    plt.grid()
    filename = 'bvp_04_solution_{:03d}.png'.format(e_num)
    plt.savefig(filename)
    print('  Graphics saved as "%s"' % (filename))
    plt.close()
#
#  Terminate.
#
    return


def bvp_04_test():

    # *****************************************************************************80
    #
    # bvp_04_test tests bvp_04.
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
    print('bvp_04_test:')
    print('  FENICS/Python version')
    print('  Solve  -u" + u = x, 0 < x < 1')
    print('  u(0) = 0, u(1) = 0 ')
    print('  Exact solution: u(x) = x - sinh(x)/sinh(1)')
    print('  Set Dirichlet boundary conditions explicitly.')

    bvp_04(e_num=8)
    bvp_04(e_num=20)
    bvp_04(e_num=50)

    print('')
    print('bvp_04_test:')
    print('  Normal end of execution.')
    print('')
    print(time.ctime(time.time()))
    return


if (__name__ == '__main__'):

    bvp_04_test()
