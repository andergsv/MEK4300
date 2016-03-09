from dolfin import *
#set_log_level(WARNING)
#set_log_active(False)
L = 4

mesh = IntervalMesh(100, 0, L)
V = FunctionSpace(mesh, 'CG', 4)
W = V * V
vh, vf = TestFunctions(W)




def start(x, on_boundary):
	return near(x[0], 0) and on_boundary

def end(x, on_boundary):
	return near(x[0], L) and on_boundary

bc = [DirichletBC(W.sub(0), Constant(0), start),
		DirichletBC(W.sub(1), Constant(0), start),
		DirichletBC(W.sub(0), Constant(1), end)]

def Newton():
	hf = Function(W)
	h, f = split(hf)

	lh = -inner(grad(h), grad(vh))*dx +2*f*h.dx(0)*vh*dx + vh*dx -h*h*vh*dx
	rh = h*vf*dx - f.dx(0)*vf*dx

	kl = lh+rh
	solve(kl == 0, hf, bc) 
	
	wiz = plot(f, interactive=False)
	wiz.write_png('Newton150')
	
def Picard():
	hf = TrialFunction(W)
	h, f = split(hf)

	hf_k = Function(W)
	hf_k = interpolate(Expression(('0', 'x[0]')), W)
	h_k, f_k = split(hf_k)

	lh = -inner( grad(h) , grad(vh) ) * dx + 2 * f_k * h.dx(0) *vh * dx + vh * dx - h_k * h * vh*dx
	rh = h*vf*dx - f.dx(0)*vf*dx

	F = lh+rh

	

	error = 1
	tol = 1.0E-12
	iter = 0
	maxiter = 200

	hf_ = Function(W)
	h_, f_ = split(hf_)
	while error>tol and maxiter>iter:		
		iter += 1
		solve(lhs(F)==rhs(F), hf_, bc)
		#h_, f_ = split(hf_)
		error = errornorm(hf_, hf_k)
		#wiz = plot(f_, interactive=False)
		#wiz.write_png('dc')
		print 'iter=%d: norm=%g' % (iter, error)
		hf_k.assign(hf_)
		h_k, f_k = split(hf_k)

	wiz = plot(f_, interactive=False)
	wiz.write_png('Picard150')
	


Newton()
#Picard()
