#! /usr/bin/env python
#
from fenics import *


def poisson_source_symbolic():

    #*****************************************************************************80
    #
    ## poisson_source_symbolic generates the right hand side of a Poisson problem.
    #
    #  Discussion:
    #
    #    We are interested in solving an equation of the form:
    #
    #      - div kappa(u,x,y) grad u = f in Omega
    #                              u = g on dOmega
    #
    #    We have chosen:
    #      u(x,y) = 2 * ( 1 + y ) / ( (3+x)^2 +(1+y)^2 )
    #      kappa(u,x,y) = 1
    #
    #    The corresponding f(x,y) is now determined.  But instead of computing it
    #    by hand, we ask the Python symbolic package to discover the formula:
    #      f(x,y) = 0
    #
    #  Licensing:
    #
    #    This code is distributed under the GNU LGPL license.
    #
    #  Modified:
    #
    #    30 October 2018
    #
    #  Author:
    #
    #    John Burkardt
    #
    import matplotlib.pyplot as plt
    import sympy as sym

    def kappa(u):
        value = 1
        return value

    x, y = sym.symbols('x[0], x[1]')
    u = 2.0 * (1.0 + y) / ((3.0 + x)**2 + (1.0 + y)**2)

    f = - sym.diff ( kappa ( u ) * sym.diff ( u, x ), x ) \
        - sym.diff ( kappa ( u ) * sym.diff ( u, y ), y )

    f = sym.simplify(f)

    f_code = sym.printing.ccode(f)
    print('')
    print('  Formula for F(x,y) computed by SymPy:')
    print('')
    print('  f(x,y) = ', f_code)
    #
    #  Terminate.
    #
    return


def poisson_source_symbolic_test():

    #*****************************************************************************80
    #
    ## poisson_source_symbolic_test tests poisson_source_symbolic.
    #
    #  Modified:
    #
    #    25 October 2018
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
    print('poisson_source_symbolic_test:')
    print('  FENICS/Python version')
    print('  Generate the right hand side function f(x,y)')
    print('  for a Poisson problem')
    print('  - div kappa(u,x,y) grad u = f(x,y) in Omega')
    print('                          u = g(x,y) on dOmega')
    print('  where we have already chosen u(x,y) and kappa(x,y).')
    print('  and want f(x,y) computed symbolically by SymPy.')

    poisson_source_symbolic()
    #
    #  Terminate.
    #
    print('')
    print('poisson_source_symbolic_test:')
    print('  Normal end of execution.')
    print('')
    print(time.ctime(time.time()))
    return


if (__name__ == '__main__'):

    poisson_source_symbolic_test()
