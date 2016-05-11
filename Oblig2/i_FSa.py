from dolfin import *
import numpy as np
import matplotlib.pyplot as plt 
parameters['reorder_dofs_serial']=False

#set_log_active(False)
beta = [1.0, 0.3, 0, -0.1, -0.18, -0.198838]

L = 6
mesh = IntervalMesh(1001, 0, L)
V = FunctionSpace(mesh, 'CG', 1)
VV = V*V
vf, vh = TestFunctions(VV)

def start(x,on_boundary):
	return near(x[0], 0) and on_boundary

def end(x,on_boundary):
	return near(x[0], L) and on_boundary

bcs = [DirichletBC(VV.sub(0),Constant(0), start), 
  		DirichletBC(VV.sub(1),Constant(0), start),
  		DirichletBC(VV.sub(1),Constant(1), end)]

def Newton(beta, guess):
	#fh = Function(VV)
	

	fh = interpolate(Expression(('%s' %guess, '%s' %guess)), VV)
	f, h = split(fh)

	lh = -inner(grad(h),grad(vh))*dx + f*h.dx(0)*vh*dx + Constant(beta)*vh*dx - Constant(beta)*h*h*vh*dx
	rh = h*vf*dx - f.dx(0)*vf*dx


	kl = lh+rh
	solve(kl == 0, fh, bcs)

	f_,h_ = fh.split(deepcopy=True)
	#h_ =project(h,V)
	
	dh = project(h.dx(0),V)
	print dh.vector().array()[0]

	#H plot
	plt.figure(1)
	plt.plot(mesh.coordinates(), h_.vector().array(), label=('beta: %s, Ig: %s' %(beta,guess)))
	#plt.axis([0,6,0,1.2])
	#plt.xlabel('L')
	#plt.title('velocity profiles')
	#plt.ylabel('f\'')
	#plt.legend(loc=4)
	plt.grid('on')
	
	#H' plot
	plt.figure(2)
	plt.plot(mesh.coordinates(), dh.vector().array(), label=('beta: %s, Ig: %s' %(beta,guess)))
	#plt.axis([0,6,0,1.4])
	#plt.xlabel('L')
	#plt.title('shear-stress profiles')
	#plt.ylabel('f\'\'')
	#plt.legend(loc=1)
	plt.grid('on')

def exa():
	for b in beta:
		Newton(b, 1)

	plt.figure(1)
	plt.axis([0,6,0,1.2])
	plt.xlabel('L')
	plt.title('velocity profiles')
	plt.ylabel('f\'')
	plt.legend(loc=4)
	plt.savefig('figs/velocityprofile.png')

	plt.figure(2)
	plt.axis([0,6,0,1.4])
	plt.xlabel('L')
	plt.title('shear-stress profiles')
	plt.ylabel('f\'\'')
	plt.legend(loc=1)
	plt.savefig('figs/shearstrssprofile.png')
	plt.show()


def exb():
	Newton(-0.1, 1)
	Newton(-0.1, 0)

	plt.figure(1)
	#plt.axis([0,6,0,1.2])
	plt.xlabel('L')
	plt.title('Two solutions for velocity profiles')
	plt.ylabel('f\'')
	plt.legend(loc=4)
	plt.savefig('figs/Twosoluvp.png')

	plt.figure(2)
	#plt.axis([0,6,0,1.4])
	plt.xlabel('L')
	plt.title('Two solutions for shear-stress profiles')
	plt.ylabel('f\'\'')
	plt.legend(loc=4)
	plt.savefig('figs/Twosolusp.png')
	plt.show()

#exa()
#exb()

Newton(-0.19884, 1)
plt.show()