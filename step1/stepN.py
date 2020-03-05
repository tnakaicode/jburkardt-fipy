#! /usr/bin/env python
#
from fenics import *

def step1 ( ):

#*****************************************************************************80
#
## step1 solves a Poisson equation with coefficient K(X).
#
#  Discussion:
#
#    STEP1:
#      Use standard finite element approach to solve a Poisson equation
#      with a coefficient function K(X).
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
#  Mesh the unit square with an 8x8 grid of triangles.
#
  m = 8
  msh = UnitSquareMesh ( m, m )
#
#  Exact solution U:
#
  U = Expression ( 'x[0]*(1-x[0])*x[1]*(1-x[1])', degree = 5 )
#
#  BLAH function K:
#
  K = Expression ( "", degree = 5 )
#
#  Exact solution right hand side F:
#
  F = Expression ( '2*x[1]*(1-x[1])+2*x[0]*(1-x[0])', degree = 5 )
#
#  Set spaces:
#    Es = Error estimator, 
#    Us = primal variable, 
#    Qs = interfacial flux.
#
  Es = FunctionSpace ( msh, 'DG', degree + 2 )
  Us = FunctionSpace ( msh, 'CG', degree + 1 )
  Qs = FunctionSpace ( msh, 'BDM', degree, restriction = 'facet' )
#
#  Set elements, needed for the MixedElement command:
#    Ee = Error estimator, 
#    Ue = primal variable, 
#    Qe = interfacial flux.
#
  Ee = FiniteElement ( 'DG',  msh.ufl_cell(), degree + 2 )
  Ue = FiniteElement ( 'CG',  msh.ufl_cell(), degree + 1 )
  Qe = FiniteElement ( 'BDM', msh.ufl_cell(), degree )
#
#  Define the mixed element, and the corresponding function space.
#
  EUQe = MixedElement ( Ee, Ue, Qe )
  EUQs = FunctionSpace ( msh, EUQe )
#
#  Get trial and test functions.
#
  ( e, u, q ) = TrialFunctions ( EUQs )
  ( y, z, r ) = TestFunctions ( EUQs )
#
#  Compute the normal vectors.
#
  n = FacetNormal ( msh )
#
#  Define the Y-inner product.
#
  yip = dot ( grad ( e ), grad ( y ) ) * dx + e * y * dx
#
#  Set up the saddle point problem:
#
#  b( (u,q), y ) = ( grad u, grad y ) - < q.n, y >
#
  b1  = dot ( grad ( u ), grad ( y ) ) * dx \
      - dot ( q ( '+' ), n ( '+' ) ) * ( y ( '+' ) - y ( '-' ) ) * dS \
      - dot ( q, n ) * y * ds
#
#  b( (z,r), e)
#
  b2 =  dot ( grad ( e ), grad ( z ) ) * dx \
      - ( e ( '+' ) - e ( '-' ) ) * dot ( r ( '+' ), n ( '+' ) ) * dS \
      - e * dot ( r, n ) * ds
#
#  Mixed form of DPG and rhs.
#
  a = yip + b1 + b2 
  b = F * y * dx
#
#  Dirichlet boundary condition.
#  DOES .sub(1) mean condition applies to second component?
#
  bc = DirichletBC ( EUQs.sub(1), U, DomainBoundary() )
#
#  Solve.
#
  euq = Function ( EUQs ) 
  solve ( a == b, euq, bcs = bc )
  e, u, q = euq.split ( deepcopy = True )
#
#  Compute errors.
#
  fmsh = refine ( refine ( msh ) )
  er = errornorm ( U, u, norm_type = 'H1', degree_rise = 3, mesh = fmsh )

  print ( '' )
  print ( '  Case of %d x %d mesh with degree %d:' % ( m, m, degree ) )
  print ( '  H1-norm of  (u - uh)   = %15.12f' % er )
  print ( '  Error estimator norm   = %15.12f' % \
    sqrt ( assemble ( ( dot ( grad ( e ), grad ( e ) ) + e * e ) * dx ) ) )
#
#  Compute the H1 projection of U (a standard Galerkin solution)
#
  uu = TrialFunction ( Us )
  vv = TestFunction ( Us )
  aa = ( dot ( grad ( uu ), grad ( vv ) ) + uu * vv ) * dx
  bb = ( F + U ) * vv * dx
  bc = DirichletBC ( Us, U, DomainBoundary() )
  pu = Function ( Us )
  solve ( aa == bb, pu, bcs = bc )
  erp = errornorm ( Us, pu, norm_type = 'H1', degree_rise = 3, mesh = fmsh )
  print ( '  Error in H1 Projection = %15.12f' % ( erp ) )
  if ( 1.0e-12 < abs ( erp ) ):
    print ( '  Error:Projection ratio = %15.12f' % ( er / erp ) )
  print ( '  H1-norm of (uh - proj) = %15.12f' % \
    errornorm ( u, pu, norm_type = 'H1', degree_rise = 0 ) )
#
#  Terminate.
#
  return

def step1_test ( ):

#*****************************************************************************80
#
## step1_test tests step1.
#
#  Modified:
#
#    23 October 2018
#
#  Author:
#
#    John Burkardt
#
  print ( '' )
  print ( 'step1:' )
  print ( '  FENICS/Python version' )
  print ( '  Discontinuous Petrov Galerkin method for the Poisson problem' )
  print ( '  on a unit square with 0 boundary conditions.' )
#
#  Terminate.
#
  print ( '' )
  print ( 'step1:' );
  print ( '  Normal end of execution.' )

if ( __name__ == '__main__' ):

  step1_test ( )
