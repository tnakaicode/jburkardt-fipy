#! /usr/bin/env python
#
from fenics import *
from mshr import *

def dpg_laplace_adapt ( ):

#*****************************************************************************80
#
## dpg_laplace_adapt applies DPG methods to the Poisson problem on the unit square.
#
#  Modified:
#
#    18 October 2018
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
""" ----------------------------------------------------------

DPG methods come with a built-in error estimator. 

This FEniCS implementation shows how to use that estimator in an adaptive
algorithm to solve

    Delta u = f    on UnitSquare
          u = 0    on boundary,

where f is peaked at a boundary point (so that near it,  we
expect the exact solution to require finer mesh points
for its resolution).

This file is part of a few graduate lectures introducing DPG
methods. It is heavily based on the undocumented FEniCS demo at

share/dolfin/demo/undocumented/adaptive-poisson/python/demo_adaptive-poisson.py

which shows a way to implement adaptivity for the standard FEM.
We have substituted the FEM by DPG, and omitted the original demo's
error estimator calculation (since one of the DPG solution components
is an error estimator).

[Disclaimer: This file worked as of May 2013 with FEniCS version
1.2.0, but it may not work on past or future versions!]

----------------------------------------Jay Gopalakrishnan  """

  import matplotlib.pyplot as plt
#
#  TOL: error tolerance for stopping adaptive iterations.
#
  TOL = 1e-4
#
#  Refine 50 % of the cells in each iteration
#        
  REFINE_RATIO = 0.50  
#
#  Maximum number of iterations.
#
  MAX_ITER = 20
#
#  Create the initial mesh.
#
  msh = UnitSquareMesh ( 2, 2 )
#
#  Set the right hand side function.
#
  source_str = "exp(-100.0*(pow(x[0], 2) + pow(x[1], 2)))"
#
#  Adaptive algorithm.
#
  for level in range ( MAX_ITER ):
#
#  Define the lowest order DPG method.
#    ER: DPG error estimator
#    CG: linear solution U
#    Q:  constant fluxes.
#
    ER = FunctionSpace ( msh, "DG", 2 )
    CG = FunctionSpace ( msh, "CG", 1 )
    Q  = FunctionSpace ( msh, "RT", 1, restriction = "facet" )
#
#  Set elements, needed for the MixedElement command:
#
    Ee = FiniteElement ( "DG", msh.ufl_cell(), 2 )
    CGe = FiniteElement ( "CG", msh.ufl_cell(), 1 )
    Qe = FiniteElement ( "RT", msh.ufl_cell(), 1 )
#
#  MixedFunctionSpace() command not available, try MixedElement instead.
#
    e3 = MixedElement ( Ee, CGe, Qe )
    X = FunctionSpace ( msh, e3 )
#
#  Define trial and test functions.
#
    ( e, u, q ) = TrialFunctions ( X )
    ( y, z, r ) = TestFunctions ( X )
#
#  Compute the normal vectors.
#
    n = FacetNormal ( msh )
#
#  Define the Y-inner product.
#
    yip = dot ( grad ( e ), grad ( y ) ) * dx + e * y * dx
#
#  Set up the saddle point problem.
#
#  b( (u,q), y ) = (grad u, grad y) - <q.n, y>
#
    b1  = dot( grad(u),grad(y) ) * dx               \
        - dot(q('+'),n('+')) * (y('+')-y('-')) * dS \
        - dot(q,n) * y * ds
#
#  b( (z,r), e)
#
    b2 =  dot( grad(e),grad(z) ) * dx \
        - (e('+')-e('-')) * dot(r('+'),n('+')) * dS \
        - e * dot(r,n) * ds
#
#  Mixed form of DPG and RHS.
#
    a = yip + b1 + b2   
    f = Expression ( source_str, degree = 10 )
    b = f * y * dx
#
#  Dirichlet conditions on the boundary.
#
    bc = DirichletBC ( X.sub(1), Constant(0.0), DomainBoundary() )      
#
#  Solve.
#
    x = Function ( X )
    solve ( a==b, x, bcs=bc )
    e, u, q = x.split()
#
#  Compute the error indicators from the DPG mixed variable e.
#
    PC = FunctionSpace ( msh, "DG", 0 )
#
#  c is a piecewise constant function.
# 
    c  = TestFunction ( PC )
    g  = assemble ( ( dot ( grad ( e ), grad ( e ) ) + e * e ) * c * dx )
#
#  Element-wise norms of e
#
    g  = g.get_local()
    E  = sqrt ( assemble ( ( dot ( grad ( e ), grad ( e ) ) + e * e ) * dx ) )
    print ( "Level %d: Estimated error E = %g (TOL = %g)" % ( level, E, TOL ) )
#
#  Check convergence.
#
    if E < TOL:
      print ( "Success, solution converged after %d iterations" % ( level ) )
      break
#
#  Mark cells for refinement.
#
    cell_markers = MeshFunction ( "bool", msh, msh.topology().dim() )
    g0 = sorted ( g, reverse = True )[int ( len ( g ) * REFINE_RATIO ) ]
    for c in cells(msh):
      cell_markers[c] = g[c.index()] > g0
#
#  Refine mesh.
#
    msh = refine ( msh, cell_markers )
#
#  Plot the mesh
#
    title_string = "Mesh for level %d" % ( level )
    plot ( msh, title = title_string )
    filename = 'mesh%d.png' % ( level )
    plt.savefig ( filename )
    print ( "Graphics: %s" % ( filename ) )
    plt.close ( )
#
#  Plot contours of the solution.
#
    title_string = "U for level %d" % ( level )
    plot ( u, mode = 'contour', title = title_string )
    filename = 'u%d.png' % ( level )
    plt.savefig ( filename )
    print ( "Graphics: %s" % ( filename ) )
    plt.close ( )
#
#  Terminate.
#
  return

def dpg_laplace_adapt_test ( ):

#*****************************************************************************80
#
## dpg_laplace_adapt_test tests dpg_laplace_adapt.
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

  print ( time.ctime ( time.time() ) )
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  print ( '' )
  print ( 'dpg_laplace_adapt_test:' )
  print ( '  FENICS/Python version' )
  print ( '  Adaptive refinement procedure applied to' )
  print ( "  Discontinuous Petrov Galerkin method for the Poisson problem" )
  print ( "  on a unit square with 0 boundary conditions." )

  dpg_laplace_adapt ( )

  print ( '' )
  print ( 'dpg_laplace_adapt_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  dpg_laplace_adapt_test ( )

