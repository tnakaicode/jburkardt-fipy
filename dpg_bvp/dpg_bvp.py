#! /usr/bin/env python
#
from fenics import *

def dpg_bvp ( ):

#*****************************************************************************80
#
## dpg_bvp applies DPG methods to a BVP on an interval.
#
#  Discussion:
#
#    This modified version of 1Dpgode.py ran on 18 October with the 
#    current version of the ever-changing FENICS.
#
#  Modified:
#
#    17 October 2018
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
#    Instead of ad-hoc technique A, one should be able to mark parts
#    of the boundary and integrate over those parts.  FEniCS has this
#    documented, but it currently seems to fail in 1D.
#
#    Instead of ad-hoc technique B, one should be able to set
#      "parameters.form_compiler.quadrature_degree"
#    to achieve high order integration.  But this facility seems to be
#    currently undocumented, so I have tried to avoid its use.
#
#    (Jay Gopalakrishnan)
#
""" ----------------------------------------------------------

This is a FEniCS implementation of a one-element (spectral)
Petrov-Galerkin method for

  u'   = f    on (0,1)
  u(0) = 0.

It was used in a graduate course at Portland to illustrate
the basic ideas of Petrov-Galerkin methods -  see the
referring example in my DPG lecture notes for details
of the method. 

We set f such that the exact solution has a layer near 1, namely

  u = ( exp(M*(x-1)) - exp(-M) ) / ( 1-exp(-M) ) 

for a large enough number M.  

[Disclaimer: This file worked as of
May 2013 with FEniCS version 1.2.0, but there is no guarantee
it will work in future versions!]

----------------------------------------Jay Gopalakrishnan  """
#
#  Make a one-element mesh.
#
  n = 1
  msh = UnitIntervalMesh ( n )
#
#  Set the maximum degree.
#
  maxp = 10
#
#  Initialize arrays to hold PG solution, L^2 project, least squares solution.
#
  sols = []
  prjs = []
  lsol = []
#
#  Set the exact solution u and corresponding right hand side f for M=30.
#
  f = Expression ( "30*exp(30*(x[0]-1))/(1-exp(-30))", degree = 30 )
  uu = Expression ( "(exp(30*(x[0]-1)) - exp(-30))/(1-exp(-30))", degree = 30 )
#
#  Set NP, the number of points where the solution will be sampled and output.
#
  np = n * ( maxp + 1 ) * 10                      
  plotmsh = UnitIntervalMesh ( np )
  ux = interpolate ( uu, FunctionSpace ( plotmsh, 'CG', 1 ) )
#
#  ad-hoc technique A to tell FEniCS to integrate the layer with precision.
#
  uhp = interpolate ( uu, FunctionSpace ( plotmsh, 'CG', maxp + 20 ) )
  fhp = interpolate ( f, FunctionSpace ( plotmsh, 'CG', maxp + 20 ) )
#
#  ad-hoc technique B to tell FEniCS to exclude left boundary point.
#
  rt = Expression ( '(x[0]<0.5 ? 0 : 1 )', degree = 10 )
#
#  Let the degree P range up to MAXP-1.
#
  for p in range ( 1, maxp ):
#
#  Compute the ideal Petrov-Galerkin solution.
#  U is the interior variable, Uhat the boundary flux variable.
#
    U = FunctionSpace ( msh, "DG", p )
    Uhat = FunctionSpace ( msh, "R", 0 )
#
#  The "MixedFunctionSpace" command has disappeared, so we have to
#  do this instead.
#
    Ue = FiniteElement ( "DG", msh.ufl_cell ( ), p );
    Uhate = FiniteElement ( "R", msh.ufl_cell ( ), 0 )
    UUhat = Ue * Uhate
#
#  X = trial space.
#
    X = FunctionSpace ( msh, UUhat )
#
#  Y = optimal test space.
#
    Y = FunctionSpace ( msh, "DG", p + 1 )

    ( u, uhat ) = TrialFunctions ( X )
    v = TestFunction ( Y )
#
#  Exclude left end point.  See "A".
#
    a = uhat * v * rt * ds - u * Dx ( v, 0 ) * dx
#
#  Error: FENICS refused "b = fhp * v * dx"
#  Precision integration: see B
#   b = fhp * v * dx
#
    b = f * v * dx
    uh = Function ( X )
    solve ( a == b, uh )
    ( u, uhat ) = uh.split ( deepcopy = True )
    ui = interpolate ( u, FunctionSpace ( plotmsh, 'CG', 1 ) )
    sols.append ( ui.vector().get_local() )
#
#  Compute the L^2 projection.
#
    u2 = TrialFunction ( U )
    v2  = TestFunction ( U )
    a2 = u2 * v2 * dx
    f2 = uu * v2 * dx
    pp = Function ( U )
    solve ( a2 == f2, pp )
    i2 = interpolate ( pp, FunctionSpace ( plotmsh,'CG',1 ) )
    prjs.append ( i2.vector().get_local() )
#
#  Compute the L^2 least squares solution.
#
    V = FunctionSpace ( msh, "CG", p ) 
    u3 = TrialFunction ( V )
    v3 = TestFunction ( V )
    a3 = Dx ( u3, 0 ) * Dx ( v3, 0 ) * dx
    f3 = f * Dx ( v3, 0 ) * dx
    def bdry ( x, on_boundary ):
      return on_boundary
#
#  Need essential boundary conditions.
#
    bc = DirichletBC ( V, rt, bdry )
    s3 = Function ( V )
    solve ( a3 == f3, s3 , bcs = bc )
    i3 = interpolate ( s3, FunctionSpace ( plotmsh, 'CG', 1 ) )
    lsol.append ( i3.vector().get_local( ) ) 
#
#  Report.
#
    print ( "  Computations with degree %d done:" % ( p ) )
    e = errornorm ( pp, u, norm_type = 'L2', degree_rise = 0, mesh = msh )
    print ( "   L^2 distance b/w projection and PG solution =%f" % ( e ) )
    e = errornorm ( uu, u, norm_type = 'L2', degree_rise = 3, mesh = plotmsh )
    print ( "   L^2 distance b/w PG soln and exact soln     =%f" % ( e ) )
    e = errornorm ( uu, s3, norm_type = 'L2', degree_rise = 3, mesh = plotmsh )
    print ( "   L^2 distance b/w least-sqr and exact soln   =%f" % ( e ) )

  print ( "Computations complete." )
#
#  Display the solutions.
#
  r = input ( "Do you want to visualize some solutions? [1/0]" )

  import numpy

  if r:
    import pylab
    x  = numpy.arange ( 0, 1+1.0/np, 1.0/np )
    Ux = ux.vector().get_local() 
    pg1 = numpy.array ( sols[0] )
    pg2 = numpy.array ( sols[int(maxp/2)] )
    pg3 = numpy.array ( sols[maxp-2] )

    figprops = dict ( figsize = ( 15, 5 ) )
    fig = pylab.figure(**figprops)

    pg = fig.add_subplot ( 1, 3, 1 )
    pg.plot ( x, pg1, 'g-.', label = 'p=1' )
    pg.plot ( x, pg2, 'b:',  label = 'p='+str(int(maxp/2)) )
    pg.plot ( x, pg3, 'r--', label = 'p='+str(maxp-1) )
    pg.plot ( x, Ux,  'k-',  label = 'Exact solution' )
    legend = pg.legend ( loc = 'upper right' )
    pylab.xlabel ( 'x' )
    pg.set_ylabel ( 'solution' )
    pg.set_title ( 'Spectral Petrov-Galerkin solutions' )

    ls = fig.add_subplot ( 1, 3, 2, sharey = pg )
    ls1 = numpy.array ( lsol[0] )      
    ls2 = numpy.array ( lsol[int(maxp/2)] ) 
    ls3 = numpy.array ( lsol[maxp-2] )
    ls.plot ( x, ls1, 'g-.', label = 'p=1' )
    ls.plot ( x, ls2, 'b:',  label = 'p='+str(int(maxp/2)) )
    ls.plot ( x, ls3, 'r--', label = 'p='+str(maxp-1) )
    ls.plot ( x, Ux,  'k-',  label = 'Exact solution' )
    legend = ls.legend ( loc = 'upper right' )
    pylab.xlabel ( 'x' )
    ls.set_title ( '$L^2$ Least-Squares solutions' )
#
#  Least squares uses p+1 for u, so subtract one degree...
#
    pr = fig.add_subplot ( 1, 3, 3, sharey = pg )
    pr1 = numpy.array ( prjs[0] )
    pr2 = numpy.array ( prjs[int(maxp/2)] )
    pr3 = numpy.array ( prjs[maxp-2] )

    pr.plot ( x, pr1, 'g-.', label = 'p=1' )
    pr.plot ( x, pr2, 'b:',  label = 'p='+str(int(maxp/2)) )
    pr.plot ( x, pr3, 'r--', label = 'p='+str(maxp-1) )
    pr.plot ( x, Ux,  'k-',  label = 'Exact solution' )
    legend = pr.legend ( loc = 'upper right' )
    pylab.xlabel ( 'x' )
    pr.set_title ( '$L^2$ projections of exact solution' )

    figfile = 'dpg_bvp.pdf'
    pylab.savefig ( figfile )
    print ( "   Saved figure to %s" % ( figfile ) )
    pylab.show()
#
#  Save data for MATLAB.
#
  r = input ( "Do you want to save solutions to Matlab? [1/0]" )

  if r:
    import scipy.io
    x  = numpy.arange ( 0, 1+1.0/np, 1.0/np )
    Ux = ux.vector().get_local() 
    pg = numpy.matrix(sols).transpose()
    pr = numpy.matrix(prjs).transpose()
    ls = numpy.matrix(lsol).transpose()
    matfile = "dpg_bvp.mat"
    scipy.io.savemat ( matfile, \
      { "PGsol":pg, "LSsol":ls, "Exact":Ux, "Proj":pr, "x":x }, \
      oned_as = 'row' )

    print ( "   Saved solution values to %s" % ( matfile ) )
    print ( "   with variables: PGsol (PG solution values)" )
    print ( "                   LSsol (LS solution values)" )
    print ( "                   Proj  (L^2 projection's values)" )
    print ( "                   Exact (exact solution values)" )
    print ( "                   x     (abscissae)." )
#
#  Terminate.
#
  return

def dpg_bvp_test ( ):

#*****************************************************************************80
#
## dpg_bvp_test tests dpg_bvp.
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

  print ( time.ctime ( time.time() ) )
#
#  Report level = only warnings or higher.
#
  level = 30
  set_log_level ( level )

  print ( '' )
  print ( 'dpg_bvp_test:' )
  print ( '  FENICS/Python version' )
  print ( '  Discontinuous Petrov Galerkin method for a boundary' )
  print ( '  value problem on the unit interval.' )

  dpg_bvp ( )

  print ( '' )
  print ( 'dpg_bvp_test:' )
  print ( '  Normal end of execution.' )
  print ( '' )
  print ( time.ctime ( time.time() ) )
  return

if ( __name__ == '__main__' ):

  dpg_bvp_test ( )
