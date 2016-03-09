from dolfin import *
from mshr import *
from numpy import log
set_log_active(False)

a = 1
b = 0.5
dpdx = -100
mu = 1
u0 = Constant(0)
meshsize=[20, 40, 80]
P = [1,2,3]

def Q_e():
	from numpy import pi 
	return pi/(4*mu)*(-dpdx)*a**3*b**3/(a**2+b**2)

er = []

for n in meshsize: 
	ellipse =  Ellipse(Point(0,0),a, b, n)
	mesh = generate_mesh(ellipse, n)

	print 'meshsize:', mesh.hmin()
	for i in range(1,4):
		V = FunctionSpace(mesh, 'CG', i)
		u = TrialFunction(V)
		v = TestFunction(V)

		def no_slip(x,on_boundary):
			return on_boundary

		bc = DirichletBC(V,u0,no_slip)

		F = -mu*inner(grad(u),grad(v))*dx == dpdx * v * dx 

		u_ = Function(V)
		solve(F,u_,bc)


		#plot(mesh, interactive=True)


		y= Expression('x[0]')
		z = Expression('x[1]')

		ue = Expression('1/(2*mu)*(-dpdx)*(a*a*b*b)/(a*a+b*b)*(1-y*y/(a*a)-z*z/(b*b))', y=y,z=z,\
		    a=Constant(a), b=Constant(b), mu=Constant(mu), dpdx=Constant(dpdx))

		uexact = project(ue, V)


		flux = u_*dx
		tot_flux = assemble(flux)



		error = errornorm(u_, uexact)
		print 'P %s :' %i, 'error:', error 
		print '      Qdif:', abs(Q_e() - tot_flux)

		er.append(error)

rp1 = [log(er[3*j]/er[3*j+3])/log(meshsize[1+j]/meshsize[j]) for j in range(2)]
rp2 = [log(er[4*j]/er[3*j+4])/log(meshsize[1+j]/meshsize[j]) for j in range(2)]
rp3 = [log(er[5*j]/er[3*j+5])/log(meshsize[1+j]/meshsize[j]) for j in range(2)]
print rp1
print rp2
print rp3
#print log(meshsize[1]/meshsize[0])


'''
uu = File('ia.pvd')
uu<<u_
'''
'''
plot(u_, title='variational')
plot(uexact,title='exact' )
plot(u_ - uexact, title='Error')
interactive()
'''

