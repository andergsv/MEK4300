from dolfin import *
set_log_level(WARNING)

def cavity_vel(ncell):
	mu = 100

	mesh = UnitSquareMesh(ncell,ncell)
	V = VectorFunctionSpace(mesh, 'CG', 2)
	Q = FunctionSpace(mesh,'CG', 1)
	VQ = V*Q
	u, p = TrialFunctions(VQ)
	v, q = TestFunctions(VQ)
	u0 = Constant((1., 0., 0.))#Expression(('1.0','0.0', '0.0'))
	no_slip =  Constant((0.,0., 0.)) # Expression(('0.0','0.0', '0.0'))

	def upper(x,on_boundary):
		return near(x[1], 1) and on_boundary

	def bottom(x, on_boundary):
		return near(x[1], 0) and on_boundary

	def right(x, on_boundary):
		return near(x[0], 1) and on_boundary

	def left(x, on_boundary):
		return near(x[0], 0) and on_boundary


	bcs =[ DirichletBC(VQ, u0, upper) , 
		   DirichletBC(VQ, no_slip, bottom), 
		   DirichletBC(VQ, no_slip, right), 
		   DirichletBC(VQ, no_slip, left)]

	F = mu*inner(grad(v), grad(u))*dx - inner(div(v), p)*dx - inner(q, div(u))*dx

	up_ = Function(VQ) 
	solve(lhs(F)==rhs(F), up_, bcs)
	u_, p_ = split(up_)


	#Streamfunc
	R = FunctionSpace(mesh, 'CG', 2)
	phi = TestFunction(R)
	psi = TrialFunction(R)

	bcstream = DirichletBC(R, 0, 'on_boundary')

	a = inner(grad(psi), grad(phi))*dx
	L = inner(u_[1].dx(0) - u_[0].dx(1), phi)*dx

	psi_ = Function(R)
	solve(a == L, psi_, bcstream)

	
	psiplot = File('iii.pvd')
	psiplot << psi_
	'''
	wiz = plot(psi_, interactive=False)
	wiz.write_png('iii')
	
	plot(psi_, interactive=True)
	'''
	minpos = psi_.vector().array().argmin()
	xco = interpolate(Expression('x[0]'), R)
	yco = interpolate(Expression('x[1]'), R)
	print 'Meshsize:', mesh.hmin()
	print 'Vortex center at x: %s and y: %s' %(xco.vector().array()[minpos], yco.vector().array()[minpos])
	
cavity_vel(100)



