'''
U = 1
L = 2
Re = 0.01
'''
from dolfin import *
set_log_level(WARNING)


L = 2
mu = 200


def top(x, on_boundary):
	return near(x[1], 0.5*L) and on_boundary

def step(x, on_boundary):
	return near(x[0], 0.5*L) and on_boundary

def bottom(x, on_boundary):
	return (near(x[1], 0.1*L) or near(x[1], 0)) and on_boundary

mesh = Mesh('meshiv.xml')

def Stokes_vel(vel, mesh):
	V = VectorFunctionSpace(mesh, 'CG', 2)
	Q = FunctionSpace(mesh,'CG', 1)
	VQ = V*Q
	u, p = TrialFunctions(VQ)
	v, q = TestFunctions(VQ)
	utop = Constant((vel, 0.))
	no_slip = Constant((0., 0.))


	bctop = DirichletBC(VQ.sub(0), utop, top)
	bcbottom = DirichletBC(VQ.sub(0), no_slip, bottom)
	bcstep = DirichletBC(VQ.sub(0), no_slip, step)

	#bcbottom = DirichletBC(VQ.sub(0), Constant((0, 0)), 'std::abs(x[1])>1e-12')
	bcs = [bctop, bcbottom, bcstep]

	F = mu*inner(grad(v), grad(u))*dx - inner(div(v), p)*dx - inner(q, div(u))*dx

	up_ = Function(VQ) 
	solve(lhs(F)==rhs(F), up_, bcs)
	u_, p_ = split(up_)
	return u_, p_


#Streamfunc
def streamfunc(u, mesh):
	u_, p_  = Stokes_vel(1., mesh)
	R = FunctionSpace(mesh, 'CG', 2)
	psiv = TestFunction(R)
	psi = TrialFunction(R)


	n = FacetNormal(mesh) 

	grad_psi = as_vector((-u[1], u[0]))

	a = -inner(grad(psi), grad(psiv))*dx 
	L = -inner(u[1].dx(0) - u[0].dx(1), psiv)*dx - psiv*dot(grad_psi,n)*ds


	psi_ = Function(R)
	solve(a == L, psi_)#, bcstream)
	normalize(psi_.vector())
	
	return R, psi_

def vortex_center(vel, d=False):
	if d ==True:
		i = 4
		mesh1 = mesh
	else:
		i = 0
		mesh1 = Mesh('meshivgrov.xml')
	while i < 5:
		u_, p_  = Stokes_vel(vel=vel, mesh=mesh1)
		R, psi = streamfunc(u_, mesh1)
		if vel > 0:
			minpos = psi.vector().array().argmin()
			val = psi.vector().array().min()
		elif vel <= 0:
			minpos = psi.vector().array().argmax()
			val = psi.vector().array().max()
		xco = interpolate(Expression('x[0]'), R)
		yco = interpolate(Expression('x[1]'), R)
		pos= [xco.vector().array()[minpos], yco.vector().array()[minpos]]
		print 'meshsize: %0.6f, minvalue: %s' %(mesh1.hmin(), val) 
		print 'Vortex center at x: %s and y: %s' %(xco.vector().array()[minpos], yco.vector().array()[minpos])
		i += 1
		mesh1 = refine(mesh1)

	if d == True:	 	
		return pos

def contour():
	u_, p_ = Stokes_vel(1., mesh)
	R, psi = streamfunc(u_, mesh)

	psiplot = File('iv.pvd')
	psiplot << psi

def flux():
	u_, p_ = Stokes_vel(1., mesh)
	def inlet(x, on_boundary):
		return near(x[0], 0) and on_boundary

	def outlet(x, on_boundary):
		return near(x[0], L) and on_boundary
	
	Inlet = AutoSubDomain(inlet)
	Outlet = AutoSubDomain(outlet)
	mf = FacetFunction('size_t', mesh)
	mf.set_all(0)
	Inlet.mark(mf,1)
	Outlet.mark(mf, 2)
	#ds = ds[mf]
	ds = Measure('ds', subdomain_data=mf)

	n = FacetNormal(mesh)
	fluxi = assemble(dot(u_,-n)*ds(1))
	fluxo = assemble(dot(u_,n)*ds(2))
	print 'Flux inn:', fluxi
	print 'Flux out:', fluxo
	print 'Flux diff:', abs(fluxi - fluxo)

def normal_stress():
	u_, p_ = Stokes_vel(-1., mesh)
	Bottom = AutoSubDomain(bottom)
	Step = AutoSubDomain(step)
	mf = FacetFunction('size_t', mesh)
	mf.set_all(0)
	Bottom.mark(mf,3)
	Step.mark(mf,4)
	ds = Measure('ds', subdomain_data=mf)

	n = FacetNormal(mesh)

	A = -assemble(p_*ds(3) + p_*ds(4))
	print 'normal stress:', A

def d():
	pos = vortex_center(1, True)
	neg =vortex_center(-1, True)
	print pos[0]-neg[0], pos[1]-neg[1]

contour()
#vortex_center(1.)
#flux()
#d()
normal_stress()

#plot(u_)
#plot(mesh, interactive=True)
#plot(psi_)
#interactive()
