from dolfin import *
from mshr import *
from numpy import sqrt, log
set_log_active(False)

dpdx = -100
mu = 1
vertices= [Point(0,0), Point(0.5,sqrt(3)/2.), Point(-0.5,sqrt(3)/2.)]
Triangle = Polygon(vertices)

meshsize= [20,40, 80]
er =[]

def Q_e():
	a = 1
	return a**4*sqrt(3)/(320*mu)*(-dpdx)

for n in meshsize:
	mesh = generate_mesh(Triangle, n)
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


		y= Expression('x[0]')
		z = Expression('x[1]')

		a = Constant(1)


		ue = Expression('-dpdx/(2*sqrt(3)*a*mu)*(z-0.5*a*sqrt(3))*(3*y*y-z*z)',y = y, z = z, \
		 a = a, mu=Constant(mu), dpdx=Constant(dpdx))
		uexact = project(ue, V)

		error = errornorm(u_, uexact)
		print 'P %s :' %i, 'error:', error 
		
		#Flux
		flux = u_*dx
		tot_flux = assemble(flux)
		print '      Qdif:', abs(tot_flux-Q_e())	
		#Save file
		'''
		uu = File('tri.pvd')
		uu<< u_
		'''

		#plot
		'''
		plot(u_, title='variational')
		plot(uexact,title='exact' )
		plot(u_ - uexact, title='Error')
		interactive()
		'''
		er.append(error)

rp1 = [log(er[3*j]/er[3*j+3])/log(meshsize[1+j]/meshsize[j]) for j in range(2)]
rp2 = [log(er[4*j]/er[3*j+4])/log(meshsize[1+j]/meshsize[j]) for j in range(2)]
rp3 = [log(er[5*j]/er[3*j+5])/log(meshsize[1+j]/meshsize[j]) for j in range(2)]
print rp1
print rp2
print rp3