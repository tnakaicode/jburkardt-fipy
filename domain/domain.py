#! /usr/bin/env python
#
from fenics import *
from mshr import *


def domain(n=10):

    # *****************************************************************************80
    #
    # domain demonstrates creating a domain by "subtraction".
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
    #  Define a region as a circle "minus" a rectangle and a square.
    #
    c1 = Point(0.0, 0.0)
    r1 = 5.0
    sw = Point(-3.0, -3.0)
    ne = Point(-2.0, 2.0)
    c2 = Point(2.0, 2.0)
    r2 = 1.0

    domain = Circle(c1, r1, segments=12) \
        - Rectangle(sw, ne) \
        - Circle(c2, r2, segments=12)
    #
    #  Mesh the region.
    #
    mesh = generate_mesh(domain, n)
    #
    #  Plot the mesh.
    #
    plot(mesh, title='domain.py: logical combination of shapes')
    filename = 'domain_mesh_{:02d}.png'.format(n)
    plt.savefig(filename)
    print('  Graphics saved as "%s"' % (filename))
    plt.close()
    #
    #  Terminate.
    #
    return


def domain_test():

    # *****************************************************************************80
    #
    # domain_test tests domain.
    #
    #  Modified:
    #
    #    23 October 2018
    #
    #  Author:
    #
    #    John Burkardt
    #

    #
    #  Report level = only warnings or higher.
    #
    level = 30
    set_log_level(level)

    print('')
    print('domain_test:')
    print('  FENICS/Python version')
    print('  Demonstrate domain creattion by subtraction.')

    domain(10)
    domain(20)
#
#  Terminate.
#
    print('')
    print('domain_test:')
    print('  Normal end of execution.')
    return


if (__name__ == '__main__'):

    domain_test()
