#! /usr/bin/env python
#
from fenics import *

def step10 ( my_grid, my_degree, my_case, my_portion ):

#*****************************************************************************80
#
## step10 plots error behavior for adapative refinement cases.
#
#  Discussion:
#
#    Use the mixed DPG method to solve a Poisson equation:
#
#      - div kappa grad u = f       in Omega = [0,1]x[0,1]
#                       u = u_exact on dOmega
#
#    where:
#
#      kappa(x,y) = 1
#
#    Do this for a sequence of refined meshes, and at the end,
#    plot the number of elements and the estimated errors.
#
#    Our refinement criterion begins by computing "g", a vector of
#    error indicator values for each cell.  We then descending sort g.
#    Then we find the index i such that g[0:i] accounts for half of the
#    error and we refine cells 0 through i.
#
#    if my_case == 1:
#
#      u_exact(x,y) = x*(1-x)*y*(1-y)
#      f(x,y)       = 2*y*(1-y)+2*x*(1-x)
#
#    if my_case == 2:
#
#      u_exact(x,y) = exp(x^2+y^2)*sin(x+y)
#      f(x,y)       = -2*exp(x^2+y^2)*(sin(x+y)*(1+2*(x^2+y^2)+2*cos(x+y)*(x+y))
#
#    if my_case == 3:
#
#      u_exact(x,y) = 2 * (1+y) / ( (3+x)^2 + (1+y)^2 )
#      f(x,y)       = 0
#
#    if my_case == 4:
#
#      u_exact(x,y) = exp ( -x^2-y^2 )
#      f(x,y)       = 4 * ( - 1 + x^2 + y^2 ) * exp ( -x^2-y^2 )
#
#    To my surprise, it seems that the degree of the Q space and element
#    must be set to 1.  Using a value of 2 results in nan values.
#    JVB, 13 November 2018.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    13 November 2018
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer MY_GRID, the number of (pairs of) elements to use in the 
#    X and Y directions.  Number of triangular elements is 2*my_grid*my_grid.
#
#    Input, integer MY_DEGREE, specifies approximation degrees.
#    Error indicators: MY_DEGREE+2
#    Primary variable: MY_DEGREE+1
#    Fluxes:           1
#
#    Input, integer MY_CASE, chooses the exact solution.
#    1 <= MY_CASE <= 4.
#
#    Input, real MY_PORTION, the proportion of the local error estimate
#    used to decide which cells to refine.
#    0.0 < MY_PORTION <= 1.
#
  import matplotlib.pyplot as plt
  import numpy as np
#
#  Report input.
#
  print ( '' )
  print ( '  Exact case #%d' % ( my_case ) )
  print ( '  %dx%d mesh on the unit square' % ( my_grid, my_grid ) )
  print ( '  Approximation degree %d:' % ( my_degree ) )
  print ( '  Error proportion %g:' % ( my_portion ) )
  print ( '' )

  if ( my_portion <= 0.0 or 1.0 < my_portion ):
    raise Exception ( '  The value of my_portion should be between 0.0 and 1.0, your value was: {}'.format ( my_portion ) )
#
#  Internal parameters:
#    LEVEL_MAX: maximum number of iterations
#
  LEVEL_MAX = 5
#
#  Mesh the unit square.
#
  my_mesh = UnitSquareMesh ( my_grid, my_grid, 'right/left' )
#
#  Define the diffusivity function kappa(x).
#
  kappa = Constant ( 1.0 )
#
#  Define expressions for the exact solution and right hand side.
#
  if ( my_case == 1 ):
    u_exact_expr = Expression ( "x[0]*(1-x[0])*x[1]*(1-x[1])", degree = 5 )
    f_expr = Expression ( "2*x[1]*(1-x[1])+2*x[0]*(1-x[0])", degree = 5 )
  elif ( my_case == 2 ):
    u_exact_expr = Expression ( "exp(pow(x[0],2)+pow(x[1],2))*sin(x[0]+x[1])", degree = 10 )
    f_expr = Expression ( "-2*exp(pow(x[0],2)+pow(x[1],2))*(sin(x[0]+x[1])*(1+2*(pow(x[0],2)+pow(x[1],2))) + 2*cos(x[0]+x[1])*(x[0]+x[1]))", degree = 10 )
  elif ( my_case == 3 ):
    u_exact_expr = Expression ( "2*(1+x[1])/(pow(3+x[0],2)+pow(1+x[1],2))", degree = 10 )
    f_expr = Expression ( "0", degree = 0 )
  elif ( my_case == 4 ):
    u_exact_expr = Expression ( "exp(-pow(x[0],2)-pow(x[1],2))", degree = 10 )
    f_expr = Expression ( "-(-4*pow(x[0],2)-4*pow(x[1],2)+4)*exp(-pow(x[0],2)-pow(x[1],2))", degree = 10 )
  else:
    raise Exception ( '  The value of my_case should be between 1 and 4, your value was: {}'.format ( my_case ) )
#
#  Initialize Numpy arrays that will hold, for each level, the number of
#  cells (elements) and the error estimate.
#
  nvec = np.zeros ( 0 )
  evec = np.zeros ( 0 )
#
#  ADAPTIVE LOOP
#
  for level in range ( LEVEL_MAX ):
#
#  Plot the current mesh.
#
    label = ( 'step10 mesh case %d level %d' % ( my_case, level ) )
    plot ( my_mesh, title = label )
    filename = ( 'step10_mesh_case%d_level%d.png' % ( my_case, level ) )
    plt.savefig ( filename )
    print ( '  Graphics saved as "%s"' % ( filename ) )
    plt.close ( )
#
#  Set spaces:
#    Es = Error estimator, 
#    Us = primal variable, 
#    Qs = interfacial flux.
#
    Es = FunctionSpace ( my_mesh, "DG",  my_degree + 2 )
    Us = FunctionSpace ( my_mesh, "CG",  my_degree + 1 )
    Qs = FunctionSpace ( my_mesh, "BDM", 1, restriction = "facet" )
#
#  Set elements, needed for the MixedElement command:
#    Ee = Error estimator element, 
#    Ue = primal variable element, 
#    Qe = interfacial flux element.
#
    Ee = FiniteElement ( "DG",  triangle, my_degree + 2 )
    Ue = FiniteElement ( "CG",  triangle, my_degree + 1 )
    Qe = FiniteElement ( "BDM", triangle, 1 )
#
#  Define the mixed element, and the corresponding function space.
#
    EUQe = MixedElement ( Ee, Ue, Qe )
    EUQs = FunctionSpace ( my_mesh, EUQe )
#
#  Extract the individual trial and test function factors from the space.
#  Here, "e" for example is a symbol for typical trial functions from Es.
#
    ( e,  u,  q  ) = TrialFunctions ( EUQs )
    ( et, ut, qt ) = TestFunctions ( EUQs )
#
#  Compute the normal vectors.
#
    n = FacetNormal ( my_mesh )
#
#  Set the components of the saddle point problem:
#
#  Define the inner product for the error estimator space.
#
    yip = dot ( grad ( e ), grad ( et ) ) * dx + e * et * dx
#
#  Set up the saddle point problem:
#
#  b( (u,q), et ) = ( kappa * grad u, grad et ) - < q.n, et >
#  This is an equation for u and q.
#
    b1 = dot ( kappa * grad ( u ), grad ( et ) ) * dx \
       - dot ( q ( '+' ), n ( '+' ) ) * ( et ( '+' ) - et ( '-' ) ) * dS \
       - dot ( q, n ) * et * ds
#
#  b( (ut,qt), e ) = ( kappa * grad ut, grad e ) - < qt.n, e >
#  This is an equation for e.
#
    b2 = dot ( kappa * grad ( ut ), grad ( e ) ) * dx \
       - dot ( qt ( '+' ), n ( '+' ) ) * ( e ( '+' ) - e ( '-' ) ) * dS \
       - dot ( qt, n ) * e * ds
#
#  Form the saddle point problem:
#
#    yip + b1 = f * et * dx
#    b2       = 0
#
    a = yip + b1 + b2 
    b = f_expr * et * dx
#
#  Indicate the Dirichlet boundary condition for the second component (u).
#
    bc = DirichletBC ( EUQs.sub(1), u_exact_expr, DomainBoundary() )
#
#  Solve the saddle point problem.
#
    euq = Function ( EUQs ) 
    solve ( a == b, euq, bcs = bc )
#
#  Extract the three components of the solution.
#  This "e, u, q" is DIFFERENT from the e, u, q above!
#  Here, "e" is the actual finite element solution function.
#
    e, u, q = euq.split ( deepcopy = True )
#
#  Compute the error estimate.
#
    cell_num = my_mesh.num_cells ( )
    error_estimate = sqrt ( assemble ( ( dot ( grad ( e ), grad ( e ) ) + e * e ) * dx ) )
    print ( '  Level %d: %d cells, error estimate is %g' % ( level, cell_num, error_estimate ) )
#
#  Add the current number of cells and error estimate to our lists.
#
    nvec = np.append ( nvec, cell_num )
    evec = np.append ( evec, error_estimate )
#
#  If this is the last level, we're done now.
#  Otherwise, compute error indicators and refine before next level.
#
    if ( LEVEL_MAX - 1 <= level ):
      print ( '' )
      print ( '  Finished with adaptive refinement sequence.' )
      break
#
#  Hs: an piecewise constant function space for error indicators.
#
    hs_degree = 0
    Hs = FunctionSpace ( my_mesh, 'DG', hs_degree )
#
#  Define a test function.
#
    ht = TestFunction ( Hs )
#
#  Define g, which is a piecewise constant mesh function whose value in each
#  cell is ||e||^2
#
    g = assemble ( ( dot ( grad ( e ), grad ( e ) ) + e * e ) * ht * dx )
#
#  Replace the mesh function by a numpy array.
#
    g = g.get_local ( )
#
#  Sum the local error estimates.
#
    g_sum = np.sum ( g )
#
#  Make a descending sorted copy of G.
#
    g_sorted = sorted ( g, reverse = True )
#
#  Find the size of the smallest error you still have to add to make
#  my_portion * g_sum.
#
    g_part = 0.0
    g_min = np.Inf

    for i in range ( 0, len ( g ) ):
      if ( g_part < my_portion * g_sum ):
        g_min = g[i]
        g_part = g_part + g_min
      else:
        break
#
#  Create a cell marker array.
#
    cell_markers = MeshFunction ( 'bool', my_mesh, my_mesh.topology().dim() )
#
#  Mark and count cells whose error norm is at least g_min.
#
    refined_num = 0

    for c in cells ( my_mesh ):
      if ( g[c.index()] >= g_min ):
        cell_markers[c] = True
        refined_num = refined_num + 1
      else:
        cell_markers[c] = False

    print ( '  Refining %d cells out of %d total' % ( refined_num, len ( g ) ) )

    if ( refined_num == 0 ):
      print ( '' )
      print ( 'step10 - Probable error!' )
      print ( '  No cells to refine?' )
      break
#
#  Refine the mesh.
#
    my_mesh = refine ( my_mesh, cell_markers )
#
#  Convert data to logarithms base 10.
#
  nvec = np.log10 ( nvec )
  evec = np.log10 ( evec )
#
#  Plot error behavior.
#
  plt.plot ( nvec, evec, marker = "o" )
  plt.title ( 'step10 error decay case %d' % ( my_case ) )
  plt.grid ( True )
  plt.xlabel ( '<-- Log10 ( Cells ) -->' )
  plt.ylabel ( '<-- Log10 ( Error ) -->' )
  filename = ( 'step10_error_decay_case%d.png' % ( my_case ) )
  plt.savefig ( filename )
  print ( '  Graphics saved as "%s"' % ( filename ) )
  plt.close ( )
#
#  Terminate.
#
  return

def step10_test ( ):

#*****************************************************************************80
#
## step10_test tests step10.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    07 November 2018
#
#  Author:
#
#    John Burkardt
#
  import time

  print ( time.ctime ( time.time() ) )
  print ( '' )
  print ( 'step10_test:' )
  print ( '  FENICS/Python version' )
  print ( '  p-Laplacian investigation.' )
  print ( '  Display error decay for adaptive refinement.' )
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  for my_case in range ( 1, 5 ):
    my_grid = 4
    my_degree = 1
    my_portion = 0.5
    step10 ( my_grid, my_degree, my_case, my_portion )
#
#  Terminate.
#
  print ( '' )
  print ( 'step10_test:' );
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )

if ( __name__ == '__main__' ):

  step10_test ( )
