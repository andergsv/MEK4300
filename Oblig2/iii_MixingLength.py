from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
parameters['reorder_dofs_serial']=False

vstar = 0.05
Ret = 1000.
nu = vstar/Ret
k = 0.41 
A = 26.

mesh = IntervalMesh(int(0.5*Ret +1),0,0.5)
#a = mesh.coordinates()
#x[:,0] = np.arctan(pi*(x[:,0])) / np.arctan(pi)
a = 1
while a*vstar/nu >= 1:
	yp = a*vstar/nu
	a = a/1.001
print 0.5/a 
#plot(mesh,interactive=True)
ypluss = Expression('x[0]*vstar/nu', vstar = vstar, nu = nu)
l = Expression('k*x[0]*(1-exp(-ypluss/A))', k = k, ypluss = ypluss, A = A)

V = FunctionSpace(mesh, 'CG', 1)


u = Function(V) #interpolate(Expression('x[0]'),V)#
v = TestFunction(V) 

def zero(x, on_boundary):
	return near(x[0], 0) and on_boundary
def one(x, on_boundary):
	return near(x[0], 1) and on_boundary

u_ = Function(V)

bcs = [DirichletBC(V, Constant(0), zero)]#,DirichletBC(V, Constant(0), lower)]

F = -Constant(nu)*inner(u.dx(0),v.dx(0))*dx + 2*vstar*vstar * v*dx - l*l*inner(abs(u.dx(0))*u.dx(0), v.dx(0))*dx


solve(F==0, u, bcs)
#plot(u,interactive=True)
def a():
	plt.plot(u.vector().array(), mesh.coordinates(), label='velocity profile')
	plt.title('Velocity profile')
	plt.grid('on')
	plt.ylabel('channel height')
	plt.legend(loc='upper left')
	plt.xlabel('u(y)')
	plt.show()


def b(B):

	yp = lambda y: y*vstar/nu 
	up5 = []#np.zeros(len(mesh.coordinates()))# np.linspace(0, 0.5, len(u.vector().array()))
	up30 = []
	coord5 = []
	coord30 = []
	j = 0
	for i in mesh.coordinates():
		if yp(i) < 5.:
			up5.append(vstar*yp(i))
			coord5.append(i)
		elif yp(i) > 30:
			up30.append(vstar*(1./k*np.log(yp(i)) + B))
			coord30.append(i)
		else:
			None
		j+=1

	if B == 4.: #Quick-fix
		plt.plot(u.vector().array(), mesh.coordinates(), label=' numerical velocity profile')
		plt.plot(up5, coord5,'o', label='y+ < 5')
	plt.plot(up30, coord30, '--', label='y+ > 30, B=%s'%B)
	
	plt.title('Velocity profile')
	plt.grid('on')
	plt.ylabel('channel height')
	plt.legend(loc='upper left')
	plt.xlabel('u(y)')
	#plt.show()

	ll = lambda y: (k*y*(1-np.exp(-yp(y)/A)))**2
	dudy = project(abs(u.dx(0)), V)
	print 'nut at center:', ll(mesh.coordinates()[-1])*dudy.vector().array()[-1]
	#print ll*dudy[-1]

a()
B = [4., 5., 5.5, 6.]

for n in B:
	b(n)
plt.show()
