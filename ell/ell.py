#! /usr/bin/env python
#
from fenics import *


def ell():

    # *****************************************************************************80
    #
    # ell solves a boundary value problem on the L-shaped region.
    #
    #  Discussion:
    #
    #    - u_xx - u_yy = f(x,y) inside L-shaped region
    #    u on boundary = 0
    #
    #    Right hand side f(x,y):
    #      rhs_fn(x,y) = e^(-(x-1)^2-(y-1)^2)
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
    #    John Burkardt
    #
    import matplotlib.pyplot as plt
#
#  Define the function used for the Dirichlet boundary conditions.
#
    ubc = Expression("0.0", degree=0)
#
#  Define the right hand side of the PDE.
#
    f = Expression("exp ( - pow ( x[0] - 1.0, 2 ) - pow ( x[1] - 1.0, 2 ) )",
                   degree=10)
#
#  Define the location of Dirichlet boundary points.
#  (Here, we say, if it's on the boundary, it's a Dirichlet point.)
#

    def x_boundary(x, on_boundary):
        return on_boundary
#
#  Read the mesh information from an XML file.
#
    mesh = Mesh('ell_mesh.xml')
#
#  Plot the mesh.
#
    plot(mesh, title='L-shaped Region Mesh')
    filename = 'ell_mesh.png'
    plt.savefig(filename)
    print('  Graphics saved as "%s"' % (filename))
    plt.close()
#
#  Define the function space (piecewise linear functions)
#
    V = FunctionSpace(mesh, "Lagrange", 1)
#
#  Define the test (perturbation) and trial (solution basis) functions.
#
    psii = TestFunction(V)
    psij = TrialFunction(V)
#
#  Define the system matrix.
#
    A = inner(grad(psii), grad(psij)) * dx
#
#  Define the system right hand side.
#
    RHS = psii * f * dx
#
#  Define the boundary conditions.
#
    bc = DirichletBC(V, ubc, x_boundary)
#
#  Indicate that the solution is a function in V.
#
    u = Function(V)
#
#  Solve the system.
#
    solve(A == RHS, u, bc)
#
#  Plot the solution.
#
    plot(u, title='L-shaped Region Solution')
    filename = 'ell_solution.png'
    plt.savefig(filename)
    print('  Graphics saved as "%s"' % (filename))
    plt.close()
#
#  Write the solution to an XML file.
#
    solution_file = File("ell_solution.xml")
    solution_file << u

    return


def ell_test():

    # *****************************************************************************80
    #
    # ell_test tests ell.
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
    print('ell_test:')
    print('  FENICS/Python version')
    print('  Boundary value problem on L-shaped region.')
    print('  -del^2 u = f(x,y) in L')
    print('  u(x,y) = 0 on dL')
    print('  f(x,y) = e^(-(x-1)^2-(y-1)^2)')

    ell()
#
#  Terminate.
#
    print('')
    print('ell_test:')
    print('  Normal end of execution.')
    print('')
    print(time.ctime(time.time()))
    return


if (__name__ == '__main__'):
    
    ell_test()
