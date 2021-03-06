#! /usr/bin/env python
#
from fenics import *


def bvp_05(e_num=10):

    # *****************************************************************************80
    #
    # bvp_05 solves a boundary value problem.
    #
    #  Discussion:
    #
    #    -u'' + u = x, 0 < x < 1
    #    u(0) = 0,
    #    u'(1) = 1-cosh(1)/sinh(1)
    #
    #    Exact solution is u(x) = x - sinh(x)/sinh(1)
    #
    #    This function has the same solution as bvp_01.py and bvp_04.py;
    #    However, here a Neumann condition on the derivative is imposed at one endpoint.
    #    This affects two portions of the problem:
    #    * the definition of the boundary conditions now drops the right endpoint;
    #    * the definition of the right hand side, where a source term must be added.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    21 October 2018
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
#  Define the right hand side.
#
    f = Expression("x[0]", degree=10)
#
#  Define the exact solution.
#
    ue = Expression("x[0] - sinh ( x[0] ) / sinh ( 1.0 )", degree=10)
#
#  Define the Neumann boundary function.
#
    g = Expression("1.0 - cosh(x[0]) / sinh(1.0)", degree=10)
#
#  Function Space
#
    V = FunctionSpace(mesh, "CG", 1)
#
#  Left boundary.
#

    def left_boundary(x, on_boundary):
        return on_boundary and (abs(x[0] - 0.0) < DOLFIN_EPS)

    u_left = ue(0.0)
    bc = DirichletBC(V, u_left, left_boundary)
#
#  Define the weak/variational form
#  Notice that "ds" is used to multiply the Neumann boundary source term.
#
    psii = TrialFunction(V)
    psij = TestFunction(V)

    A = inner(grad(psii), grad(psij)) * dx + psii * psij * dx
    RHS = f * psij * dx + g * psij * ds
#
#  Solve
#
    u = Function(V)

    solve(A == RHS, u, bc)
#
#  Plot the solution.
#
    plot(u, title='bvp_05 solution')
    plt.grid()
    filename = 'bvp_05_solution_{:03d}.png'.format(e_num)
    plt.savefig(filename)
    print('  Graphics saved as "%s"' % (filename))
    plt.close()

    return


def bvp_05_test():

    # *****************************************************************************80
    #
    # bvp_05_test tests bvp_05.
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
    print('bvp_05_test:')
    print('  FENICS/Python version')
    print('  Solve  -u" + u = x, 0 < x < 1')
    print('  u(0) = 0')
    print('  u''(1) = 1-cosh(1)/sinh(1)')
    print('  Exact solution: u(x) = x - sinh(x)/sinh(1)')
    print('  Notice Neumann condition on the right.')

    bvp_05(e_num=8)
    bvp_05(e_num=10)
    bvp_05(e_num=20)
    bvp_05(e_num=50)

    print('')
    print('bvp_05_test:')
    print('  Normal end of execution.')
    print('')
    print(time.ctime(time.time()))
    return


if (__name__ == '__main__'):
    
    bvp_05_test()
