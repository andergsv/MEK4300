from dolfin import *
from mshr import *
import numpy as np
set_log_active(False)



N = 60
nc = 100


dt = 0.001
T = 8

nu = Constant(0.001)
rho = Constant(1)
#Re = Constant(20)

cyldia = 0.1
L = 2.2
H = cyldia + 0.15 + 0.16

rec = Rectangle(Point(0,0), Point(L, H))
cyl = Circle(Point(0.15+cyldia/2., 0.15+cyldia/2.), cyldia/2., nc)
def coeffs(u, p, normal, surface, Um, Dia, mf):
	rho = 1
	Ua = 2*Um/3.
	normal = -normal
	ds = Measure("ds", subdomain_id=surface,subdomain_data=mf)
	
	'''
	n1 = as_vector((1., 0))
	n2 = as_vector((0, 1.))
	nx = dot(normal ,n1)
	ny = dot(normal, n2)
	nt = as_vector((ny, -nx))
	
	Dr = (rho*nu*dot(grad(dot(u,nt)),normal)*ny - p*nx)*ds
	Li = -(rho*nu*dot(grad(dot(u,nt)),normal)*nx + p*ny)*ds


	'''
	eps = 0.50*(nabla_grad(u) + nabla_grad(u).T)
	o = 2*nu*rho*eps - p*Identity(2)

	traction = dot(o, normal)

	Dr = traction[0]*ds#(surface)
	Li = traction[1]*ds#(surface)
	
	
	Fd = assemble(Dr)
	Fl = assemble(Li)

	cd = 2.*Fd/(rho*Ua**2*Dia)
	cl = 2.*Fl/(rho*Ua**2*Dia)#(rho*Um*Um*Dia)
	
	return cd, cl

def DeltaP(up, task):
	if task==2:
		Pa = up(Point(0.15,0.2))
		Pe = up(Point(0.25,0.2))
		delP =  Pa - Pe#0.1172
		return delP
	else:
		Pa = up(Point(0.15,0.2))[2]
		Pe = up(Point(0.25,0.2))[2]
		delP =  Pa - Pe#0.1172
		return delP

mesh = generate_mesh(rec - cyl, N)
#mesh = Mesh('ii.xml')

for i in range(3):
	if i == 0:
		None
	else:
		mesh = refine(mesh)
	#plot(mesh, interactive=True)
	V = VectorFunctionSpace(mesh, 'CG', 2)
	Q = FunctionSpace(mesh, 'CG', 1)
	W = V*Q

	mf = FacetFunction('size_t', mesh)

	n = FacetNormal(mesh)

	mf.set_all(0)

	class Walls(SubDomain):
		def inside(self, x, on_boundary):
			return on_boundary and near(x[1], H) or near(x[1],0) 

	class Cylinder(SubDomain):
		def inside(self, x, on_boundary):
			return on_boundary and not near(x[1], H) and not near(x[1],0) \
							   and not near(x[0], 0) and not near(x[0],L)

	class Inlet(SubDomain):
		def inside(self, x, on_boundary):
			return on_boundary and near(x[0], 0) 

	class Outlet(SubDomain):
		def inside(self, x, on_boundary):
			return on_boundary and near(x[0], L)

	walls = Walls()
	walls.mark(mf,1)
	cylinder = Cylinder()
	cylinder.mark(mf, 2)
	inlet = Inlet()
	inlet.mark(mf,3)
	outlet = Outlet()
	outlet.mark(mf,4)
	#plot(mf, interactive=True)
	no_slip = Constant((0,0))



	Re = Constant(100)
	Um = Constant(1.5)

	bcwalls = DirichletBC(V, no_slip, walls)
	bccyl = DirichletBC(V, no_slip, cylinder)

	U = Expression(('4*Um*x[1]*(H-x[1])/(H*H)','0'),Um = Um, H = Constant(H))
	#plot(mf, interactive=True)

	bcinlet = DirichletBC(V, U, inlet)
	bcoutlet = DirichletBC(Q, Constant(0), outlet)

	bcv = [bcinlet, bccyl, bcwalls]
	bcp = bcoutlet


	u = Function(V)
	p = Function(Q)

	v = TestFunction(V)
	q = TestFunction(Q)

	#Initial Condition
	u0 = Function(V)
	p0 = Function(Q)

	u_ = Function(V)

	eps = lambda u: 0.50*(grad(u) + grad(u).T)
	o = lambda u: 2*nu*rho*eps(u) - p0*Identity(2)


	'''

	F1 = u*v*dx + dt*inner(dot(u,nabla_grad(u)), v)*dx + dt*rho*nu*inner(grad(u), grad(v))*dx + dt*inner(grad(p),v)*dx -u0*v*dx
	F2 = inner(div(u),q)*dx
	'''
	kn = dt


	F1 = rho/kn*inner((u-u0),v)*dx + rho*inner(dot(u0,nabla_grad(u0)),v)*dx  + inner(o((u+u0)/2.),eps(v))*dx \
		-nu*rho*inner(dot(grad(u0),n),v)*ds(4) + inner(p0*n, v)*ds(4) 

	F2 = kn*inner(grad(p), grad(q))*dx - kn*inner(grad(p0), grad(q))*dx + rho*inner(div(u),q)*dx

	F3 = rho*inner(u_,v)*dx - rho*inner(u,v)*dx+ kn*inner(grad(p - p0), v)*dx
	if i ==0:
		file = File('para/NSunsteady.pvd')

	t = 0
	j = 0

	outfile = open('output%s.txt'%i, 'w')
	w = False
	while t<=T:
		solve(F1 == 0, u, bcv)
		solve(F2 == 0, p, bcp)
		solve(F3 == 0, u_, bcv)


		p0.assign(p)
		u0.assign(u_)

		cd, cl = coeffs(u_, p, n, 2, 1.5, 0.1, mf)
		
		delP = DeltaP(p,2)

		#print p(Point(0.15,0.2))

		j += 1
		t += dt

		print '{:>10}%      dofs={:>10}        t={:>5}      Cd={:>10}          Cl={:>17}     delP={:>15} '.format(t/float(T)*100,W.dim(), t, cd, cl, delP)
		
		if j == 20:
			w = True
			if i == 0:
				file << u
				j = 0

		#plot(u_)#, interactive=True)
		
		if w==True:
			outfile.write('%0.7f %5.7f %5.7f \n' %(cd, cl, delP))

	outfile.close()




