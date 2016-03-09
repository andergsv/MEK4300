from dolfin import *
import numpy as np
set_log_active(False)

L = 4
mesh = IntervalMesh(50, 0, L)
V = FunctionSpace(mesh, 'CG', 1)
VV = V*V
vf, vh = TestFunctions(VV)

def start(x,on_boundary):
	return near(x[0], 0) and on_boundary

def end(x,on_boundary):
	return near(x[0], L) and on_boundary

bc = [DirichletBC(VV.sub(0),Constant(0), start), 
  		DirichletBC(VV.sub(1),Constant(0), start),
  		DirichletBC(VV.sub(1),Constant(1), end)]

def Newton():
	fh = Function(VV)
	f, h = split(fh)

	lh = -inner(grad(h),grad(vh))*dx + f*h.dx(0)*vh*dx + vh*dx - h*h*vh*dx
	rh = h*vf*dx - f.dx(0)*vf*dx


	kl = lh+rh
	solve(kl == 0, fh, bc)
	'''
	F_h = Function(VV)
	solve(lhs(kl) ==rhs(kl), F_h, bc=[bcf,bch,bce])
	'''
	wiz = plot(f, interactive=False)
	wiz.write_png("Newton148")

def Picard():
	
	fh = TrialFunction(VV)
	f, h = split(fh)	

	fh_k = Function(VV)
	f_k, h_k = split(fh_k)
	
	#h_k = interpolate(Expression('x[0]'), V)p

	lh = -inner(grad(h),grad(vh))*dx  + vh*dx - h_k*h*vh*dx + h_k.dx(0)*f*vh*dx
	rh = h*vf*dx - f.dx(0)*vf*dx

	F = lh+rh
	fh_ = Function(VV)
	
	

	error = 1.0            # error measure ||u-u_k||
	tol = 1.0E-12          # tolerance
	iter = 0               # iteration counter
	maxiter = 200           # max no of iterations allowed
	while error > tol and iter < maxiter:
	    iter += 1
	    solve(lhs(F)==rhs(F), fh_, bc) #solve(a == L, u, bcs)
	    f_, h_ = split(fh_)
	    error = errornorm(fh_, fh_k)
	    print 'iter=%d: norm=%g' % (iter, error)
	    fh_k.assign(fh_)   # update for next iteration
	    f_k, h_k = split(fh_k)

	wiz = plot(f_, interactive=False)
	wiz.write_png('Picard148')


Newton()
#Picard()