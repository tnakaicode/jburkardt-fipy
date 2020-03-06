#! /usr/bin/env python3
#
from fenics import *


def heat_steady_02():

    #*****************************************************************************80
    #
    ## heat_steady_02, 2D steady heat equation with point source.
    #
    #  Discussion:
    #
    #    Heat equation with point source:
    #
    #      - uxx - uyy = f(x,y) = 2 + 100.0 * delta(0.0,0.0)
    #      over the rectangle -1.0 <= x, y <= 1.0
    #      with zero Dirichlet boundary conditions.
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    10 January 2016
    #
    #  Author:
    #
    #    John Burkardt
    #
    import matplotlib.pyplot as plt
    #
    #  Define the mesh.
    #    The region is a rectangle with corners at ( -1,-1) and (+1,+1).
    #    The mesh uses 100 divisions in the X and Y directions.
    #
    sw = Point(-1.0, -1.0)
    ne = Point(+1.0, +1.0)
    mesh = RectangleMesh(sw, ne, 100, 100)
    #
    #  Define the function space.
    #
    V = FunctionSpace(mesh, "CG", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    #
    #  Define the bilinear form, the left hand side of the FEM system.
    #
    Auv = inner(grad(u), grad(v)) * dx
    #
    #  Define the linear form, the right hand side.
    #
    Lv = Constant(2.0) * v * dx
    #
    #  If a Dirichlet boundary condition is to be applied exactly on the
    #  boundary of the domain, we can use "DomainBoundary" as the
    #  third argument to DirichletBC.
    #
    bc = DirichletBC(V, Constant(0.0), DomainBoundary())
    #
    #  Assemble the system matrix and right hand side vector,
    #  enforcing the boundary conditions.
    #
    A, B = assemble_system(Auv, Lv, bc)
    #
    #  Modify the right hand side vector to account for a point source
    #  located at (0,0) with magnitude 100.
    #
    delta = PointSource(V, Point(0.0, 0.0), 100.0)
    delta.apply(B)
    #
    #  Solve the linear system for u.
    #
    u = Function(V)
    solve(A, u.vector(), B)
    #
    #  Plot the solution.
    #
    plot(u, title='heat_steady_02')
    filename = 'heat_steady_02.png'
    plt.savefig(filename)
    print("Saving graphics in file '%s'" % (filename))
    plt.close()
    #
    #  Terminate.
    #
    return


def heat_steady_02_test():

    #*****************************************************************************80
    #
    ## heat_steady_02_test tests heat_steady_02.
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
    print('heat_steady_02_test:')
    print('  FENICS/Python version')
    print('  Solve the heat equation over a square, with a point source.')

    heat_steady_02()
    #
    #  Terminate.
    #
    print('')
    print('heat_steady_02_test:')
    print('  Normal end of execution.')
    print('')
    print(time.ctime(time.time()))
    return


if (__name__ == '__main__'):

    heat_steady_02_test()
