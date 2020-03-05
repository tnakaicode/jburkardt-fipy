#! /usr/bin/env python
#
from fenics import *


def dpg_laplace():

    # *****************************************************************************80
    #
    # dpg_laplace applies DPG methods to the Poisson problem on the unit square.
    #
    #  Modified:
    #
    #    16 October 2018
    #
    #  Modifier:
    #
    #    John Burkardt
    #
    #  Author:
    #
    #    Jay Gopalakrishnan
    #
    #  Reference:
    #
    #    Jay Gopalakrishnan,
    #    Five lectures on DPG Methods,
    #    Spring 2013, Portland State University,
    #    arXiv:1306.0557v2 [math.NA] 28 Aug 2014.
    #
    #  Issues:
    #
    #    If you put degree=2 then the H1 best approximation error comes
    #    out LARGER than the error of a numerical solution in the same
    #    space!  This is mathematically impossible.  Is FEniCS losing some
    #    double-precision digits in some routine?  (Jay Gopalakrishnan)
    #
    """ ----------------------------------------------------------

    This is a FEniCS implementation of a DPG method for the
    Dirichlet problem

        Delta u = f    on UnitSquare
              u = g    on boundary,

    and is part of notes for graduate lectures introducing DPG methods.

    [Disclaimer: This file worked as of May 2013 with FEniCS version
    1.2.0. It may not work in past or future versions!]

    ----------------------------------------Jay Gopalakrishnan  """
    #
    #  Degree of the trial space.
    #
    degree = 1
    #
    #  Let mesh be a square, divided into an 8x8 grid of squares, triangulated.
    #
    m = 8
    msh = UnitSquareMesh(m, m)
    #
    #  Exact solution U:
    #
    U = Expression("x[0]*(1-x[0])*x[1]*(1-x[1])", degree=5)
    #
    #  Exact solution right hand side F:
    #
    F = Expression("2*x[1]*(1-x[1])+2*x[0]*(1-x[0])", degree=5)
    #
    #  Second possible exact solution U:
    #
    # U = Expression ( "exp(pow(x[0],2)+pow(x[1],2))*sin(x[0]+x[1])" )
    #
    #  Second possible right hand side F.
    #
    # F = Expression ( "-2*exp(pow(x[0],2)+pow(x[1],2))*(sin(x[0]+x[1])*(1+2*(pow(x[0],2)+pow(x[1],2))) + 2*cos(x[0]+x[1])*(x[0]+x[1]))" )
    #
    #  Set spaces:
    #    E = Error estimator,
    #    CG = primal variable,
    #    Q = interfacial flux.
    #
    E = FunctionSpace(msh, "DG", degree + 2)
    CG = FunctionSpace(msh, "CG", degree + 1)
    Q = FunctionSpace(msh, "BDM", degree, restriction="facet")
    #
    #  Set elements, needed for the MixedElement command:
    #    E = Error estimator,
    #    CG = primal variable,
    #    Q = interfacial flux.
    #
    Ee = FiniteElement("DG", msh.ufl_cell(), degree + 2)
    CGe = FiniteElement("Lagrange", msh.ufl_cell(), degree + 1)
    Qe = FiniteElement("BDM", msh.ufl_cell(), degree)
    #
    #  Latest FENICS removed MixedFunctionSpace() command, thank you very much...
    #  The work-around was not easy to find.
    #
    e3 = MixedElement(Ee, CGe, Qe)
    X = FunctionSpace(msh, e3)
    #
    #  Get trial and test functions.
    #
    (e, u, q) = TrialFunctions(X)
    (y, z, r) = TestFunctions(X)
    #
    #  Compute the normal vectors.
    #
    n = FacetNormal(msh)
    #
    #  Define the Y-inner product.
    #
    yip = dot(grad(e), grad(y)) * dx + e * y * dx
    #
    #  Set up the saddle point problem:
    #
    #  b( (u,q), y ) = ( grad u, grad y ) - < q.n, y >
    #
    b1 = dot(grad(u), grad(y)) * dx \
        - dot(q('+'), n('+')) * (y('+') - y('-')) * dS \
        - dot(q, n) * y * ds
    #
    #  b( (z,r), e)
    #
    b2 = dot(grad(e), grad(z)) * dx \
        - (e('+') - e('-')) * dot(r('+'), n('+')) * dS \
        - e * dot(r, n) * ds
    #
    #  Mixed form of DPG and rhs.
    #
    a = yip + b1 + b2
    b = F * y * dx
    #
    #  Dirichlet boundary condition.
    #
    bc = DirichletBC(X.sub(1), U, DomainBoundary())
    #
    #  Solve.
    #
    x = Function(X)
    solve(a == b, x, bcs=bc)
    e, u, q = x.split(deepcopy=True)
    #
    #  Compute errors.
    #
    fmsh = refine(refine(msh))
    er = errornorm(U, u, norm_type='H1', degree_rise=3, mesh=fmsh)

    print("")
    print("  Case of %d x %d mesh with degree %d:" % (m, m, degree))
    print("  H1-norm of  (u - uh)   = %15.12f" % er)
    print("  Error estimator norm   = %15.12f" %
          sqrt(assemble((dot(grad(e), grad(e)) + e * e) * dx)))
    #
    #  Compute the H1 projection of U (a standard Galerkin solution)
    #
    uu = TrialFunction(CG)
    vv = TestFunction(CG)
    aa = (dot(grad(uu), grad(vv)) + uu * vv) * dx
    bb = (F + U) * vv * dx
    bc = DirichletBC(CG, U, DomainBoundary())
    pu = Function(CG)
    solve(aa == bb, pu, bcs=bc)
    erp = errornorm(U, pu, norm_type='H1', degree_rise=3, mesh=fmsh)
    print("  Error in H1 Projection = %15.12f" % (erp))
    if (1.0e-12 < abs(erp)):
        print("  Error:Projection ratio = %15.12f" % (er / erp))
    print("  H1-norm of (uh - proj) = %15.12f" %
          errornorm(u, pu, norm_type='H1', degree_rise=0))
    #
    #  Terminate.
    #
    return


def dpg_laplace_test():

    # *****************************************************************************80
    #
    # dpg_laplace_test tests dpg_laplace.
    #
    #  Modified:
    #
    #    29 October 2018
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
    print('dpg_laplace_test:')
    print('  FENICS/Python version')
    print("  Discontinuous Petrov Galerkin method for the Poisson problem")
    print("  on a unit square with 0 boundary conditions.")

    dpg_laplace()

    print('')
    print('dpg_laplace_test:')
    print('  Normal end of execution.')
    print('')
    print(time.ctime(time.time()))
    return


if (__name__ == '__main__'):

    dpg_laplace_test()
