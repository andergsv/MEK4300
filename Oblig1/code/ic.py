from dolfin import *
from mshr import *
from numpy import sqrt, log, sinh, pi, exp
set_log_active(False)

dpdx = -100
mu = 1
a = 1
b = 0.2
c = 0.5
meshsize = [20,40, 80]

def Q_e():
	F = ( a**2 - b**2 + c**2) / ( 2 * c )
	M = sqrt( F**2 - a**2 )
	alpha = 0.5 * log( ( F + M ) / ( F - M ) )
	beta = 0.5 * log( ( F - c + M ) / ( F - c - M ) )
	Q = (a**4 - b**4 - 4 * c**2 * M**2 / (beta-alpha))
	n = 2
	Q1 = 8*c**2*M**2* exp(-(beta+alpha))/sinh(beta - alpha)
	err = 1
	Q -= Q1 
	while err>1e-14:
		Qsum = 8*c**2*M**2* n*exp(-n*(beta+alpha))/sinh(n*beta - n* alpha) 
		err = abs(Qsum - Q1)
		Q -= Qsum
		Q1 = Qsum
		n += 1
	return pi/(8*mu)*(-dpdx)*Q

for n in meshsize:
	o_circ = Circle(Point(0,0),a, n)
	i_circ = Circle(Point(c,0),b, n)
	mesh = generate_mesh(o_circ-i_circ, n)
	print 'meshsize:', mesh.hmin()
	for i in range(1,4):
		V = FunctionSpace(mesh, 'CG', i)
		u = TrialFunction(V)
		v = TestFunction(V)

		def no_slip(x,on_boundary):
			return on_boundary

		bc = DirichletBC(V, Constant(0),no_slip)

		F = -mu*inner(grad(u),grad(v))*dx == dpdx * v * dx 

		u_ = Function(V)

		solve(F,u_,bc)

		#Flux
		flux= u_*dx
		tot_flux = assemble(flux)
		Qdif = abs(tot_flux - Q_e())
		print 'P %s' %i, 'Qdif:', Qdif
		#print tot_flux, Q_e()
		#Save FIle
		'''
		uu= File('Anul.pvd')
		uu<<u_
		'''

		#PLot
		'''
		plot(mesh,title='mesh')
		plot(u_, title='variational')
		interactive()
		'''	