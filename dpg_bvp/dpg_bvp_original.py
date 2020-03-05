""" ----------------------------------------------------------

This is a FEniCS implementation of a one-element (spectral)
Petrov-Galerkin method for

  u'   = f    on (0,1)
  u(0) = 0.

It was used in a graduate course at Portland to illustrate
the basic ideas of Petrov-Galerkin methods -  see the
referring example in my DPG lecture notes for details
of the method. We set f such that the exact solution has a layer
near 1, namely

  u = ( exp(M*(x-1)) - exp(-M) ) / ( 1-exp(-M) ) 

for a large enough number M.  [Disclaimer: This file worked as of
May 2013 with FEniCS version 1.2.0, but there is no guarantee
it will work in future versions!]

----------------------------------------Jay Gopalakrishnan  """

from dolfin import *

# make a one-element mesh and init
n = 1
msh = UnitIntervalMesh(n)              
maxp = 10                              # max degree
sols = []                              # values of PG solution
prjs = []                              # values of L^2 projection
lsol = []                              # values of least-squares sol

# set exact solution u and corresponding rhs f for M=30
f  = Expression("30*exp(30*(x[0]-1))/(1-exp(-30))")
uu = Expression('(exp(30*(x[0]-1)) - exp(-30))/(1-exp(-30))')

# num points where solution will be sampled and output
np = n*(maxp+1)*10                      
plotmsh = UnitIntervalMesh(np)
ux = interpolate(uu,FunctionSpace(plotmsh,'CG',1))

# ad-hoc technique A to tell FEniCS to integrate the layer with precision
uhp= interpolate(uu,FunctionSpace(plotmsh,'CG',maxp+20))
fhp= interpolate(f,FunctionSpace(plotmsh,'CG',maxp+20))

# ad-hoc technique B to tell FEniCS to exclude left boundary point
rt = Expression('(x[0]<0.5 ? 0 : 1 )')

for p in range(1,maxp):                # for degrees from 1 to maxp:

    # compute ideal Petrov-Galerkin solution
    U = FunctionSpace(msh, "DG", p)    # interior variable
    Uhat = FunctionSpace(msh, "R",0)   # boundary 'flux' variable
    X = MixedFunctionSpace([U,Uhat])   # trial space
    Y = FunctionSpace(msh, "DG", p+1)  # optimal test space
    (u,uhat) = TrialFunctions(X)
    v = TestFunction(Y)
    a = uhat*v*rt*ds  - u*Dx(v,0)*dx   # exclude left end pt: see A
    b = fhp*v*dx                       # precise integration: see B
    uh = Function(X)
    solve( a==b, uh )                  # solve
    (u,uhat) = uh.split(deepcopy=True)
    ui = interpolate(u,FunctionSpace(plotmsh,'CG',1))
    sols.append( ui.vector().array() ) # store for output

    # Issues:
    #
    # Instead of ad-hoc technique A, one should be able to mark parts
    # of the boundary and integrate over those parts. FEniCS has this
    # documented, but it currently seems to fail in 1D.
    #
    # Instead of ad-hoc technique B, one should be able to set
    #     "parameters.form_compiler.quadrature_degree"
    # to achieve high order integration. But this facility seems to be
    # currently undocumented, so I have tried to avoid its use.
    
    # compute L^2 projection
    u2 = TrialFunction(U)
    v2  = TestFunction(U)
    a2 = u2*v2*dx
    f2 = uhp*v2*dx
    pp = Function(U)
    solve( a2 == f2, pp )
    i2 = interpolate(pp,FunctionSpace(plotmsh,'CG',1))
    prjs.append( i2.vector().array() ) 

    # compute L^2 least squares solution
    V = FunctionSpace(msh, "CG", p) 
    u3 = TrialFunction(V)
    v3 = TestFunction(V)
    a3 = Dx(u3,0) * Dx(v3,0) * dx
    f3 = fhp*Dx(v3,0) * dx
    def bdry(x,on_boundary):
        return on_boundary
    bc = DirichletBC(V,rt,bdry)        # need essential bc
    s3 = Function(V)
    solve( a3 == f3, s3 , bcs=bc)
    i3 = interpolate(s3,FunctionSpace(plotmsh,'CG',1))
    lsol.append( i3.vector().array() ) 
    
    # report
    print " >Computations with degree %d done:" % p
    e = errornorm(pp,u,norm_type='L2',degree_rise=0,mesh=msh)
    print "   L^2 distance b/w projection and PG solution =%f"% e
    e = errornorm(uu,u,norm_type='L2',degree_rise=3,mesh=plotmsh)
    print "   L^2 distance b/w PG soln and exact soln     =%f"% e
    e = errornorm(uu,s3,norm_type='L2',degree_rise=3,mesh=plotmsh)
    print "   L^2 distance b/w least-sqr and exact soln   =%f"% e

print "Computations complete."

# Plot, output etc. (Make sure you have installed numpy, pylab,
# and scipy if you want the rest of the code to work!)
r = input("Do you want to visualize some solutions? [1/0] ")
import numpy
if r:
    import pylab
    x  = numpy.arange(0,1+1.0/np,1.0/np)
    Ux = ux.vector().array() 
    pg1= numpy.array(sols[0])
    pg2= numpy.array(sols[maxp/2])
    pg3= numpy.array(sols[maxp-2])

    figprops = dict(figsize=(15, 5))
    fig = pylab.figure(**figprops)
    pg = fig.add_subplot(1, 3, 1)
    pg.plot(x,pg1,'g-.',x,pg2,'b:',x,pg3,'r--', x,Ux,'k-')
    pylab.legend(('p=1', 'p='+str(maxp/2), 'p='+str(maxp-1),\
                  'Exact solution'),'upper left')
    pylab.xlabel('x')
    pg.set_ylabel('solution')
    pg.set_title('Spectral Petrov-Galerkin solutions')

    ls = fig.add_subplot(1, 3, 2, sharey=pg)
    ls1= numpy.array(lsol[0])      
    ls2= numpy.array(lsol[maxp/2]) 
    ls3= numpy.array(lsol[maxp-2])
    pylab.plot(x,ls1,'g-.',x,ls2,'b:',x,ls3,'r--', x,Ux,'k-')
    pylab.legend(('p=1', 'p='+str(maxp/2), 'p='+str(maxp-1),\
                  'Exact solution'),'upper left')             
    pylab.xlabel('x')
    ls.set_title('$L^2$ Least-Squares solutions')
    
    pr = fig.add_subplot(1, 3, 3, sharey=pg)
    pr1= numpy.array(prjs[0])           # least squares uses p+1 for u,
    pr2= numpy.array(prjs[maxp/2])      # so subtract one degree
    pr3= numpy.array(prjs[maxp-2])
    pylab.plot(x,pr1,'g-.',x,pr2,'b:',x,pr3,'r--', x,Ux,'k-')
    pylab.legend(('p=1', 'p='+str(maxp/2), 'p='+str(maxp-1),\
                  'Exact solution'),'upper left')
    pylab.xlabel('x')
    pr.set_title('$L^2$ projections of exact solution')

    figfile = 'pglsprj.pdf'
    pylab.savefig(figfile)
    print "   Saved figure to %s"% figfile 
    pylab.show()

r = input("Do you want to save solutions to Matlab? [1/0] ")
if r:
    import scipy.io
    x  = numpy.arange(0,1+1.0/np,1.0/np)
    Ux = ux.vector().array() 
    pg = numpy.matrix(sols).transpose()
    pr = numpy.matrix(prjs).transpose()
    ls = numpy.matrix(lsol).transpose()
    matfile = "pglsprj.mat"
    scipy.io.savemat(matfile,{"PGsol":pg,"LSsol":ls,"Exact":Ux,\
                              "Proj":pr, "x":x},oned_as='row')
                              
    print "   Saved solution values to %s"% matfile
    print "   with variables: PGsol (PG solution values)"
    print "                   LSsol (LS solution values)"
    print "                   Proj  (L^2 projection's values)"
    print "                   Exact (exact solution values)"
    print "                   x     (abscissae)."

