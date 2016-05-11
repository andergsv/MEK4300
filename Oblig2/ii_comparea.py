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




rec = Rectangle(Point(0,0), Point(L, H))
cyl = Circle(Point(0.15+cyldia/2., 0.15+cyldia/2.), cyldia/2., nc)
#mesh = generate_mesh(rec - cyl, N)
mesh = Mesh('ii.xml')


print 'dofs      Cd       Cl       La     DeltaP'
for i in range(2):
	#plot(mesh, interactive=True)
	if i == 0:
		None
	else:
		mesh = refine(mesh)
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


	Re = Constant(20)
	Um = Constant(0.3)

	up = Function(W)
	v, q = TestFunctions(W)

	u, p = split(up)

	bcwalls = DirichletBC(W.sub(0), no_slip, walls)
	bccyl = DirichletBC(W.sub(0), no_slip, cylinder)

	U = Expression(('4*Um*x[1]*(H-x[1])/(H*H)','0'),Um = Um, H = Constant(H))
	#plot(mf, interactive=True)


	bcinlet = DirichletBC(W.sub(0), U, inlet)
	bcoutlet = DirichletBC(W.sub(1), Constant(0), outlet)

	bcs = [bcinlet, bccyl, bcwalls, bcoutlet]

	F1 = inner(dot(u,nabla_grad(u)), v)*dx + rho*nu*inner(grad(u), grad(v))*dx + inner(grad(p),v)*dx
	F2 = inner(div(u),q)*dx

	F = F1 - F2

	solve(F == 0, up, bcs)


	cd, cl = coeffs(u, p, n, 2, 0.3, 0.1, mf)



	xe = 0.25
	x = np.linspace(xe, L, 10001)

	for i in x:
		a = up(Point(i,0.2))[0]
		#b = up(Point(i,0.2))[1]
		if a > 1e-12 :
			xr = i
			break
	La = xr-xe #0.0845
	delP = DeltaP(up,1)

	print '%s   %.4f   %.4f   %.4f   %.4f' %(W.dim(), cd, cl, La, delP)
	'''
	print 'Cd:', cd,'Cl:', cl 
	print 'La:', La
	print 'delP:', delP
	print W.dim()
	'''
	#plot(p, interactive=True)
	#plot(u, interactive=True)
	#uplot = plot(u, interactive=False)
	#uplot.write_png('u')
	#pplot = plot(p, interactive=False)
	#pplot.write_png('p')







#oppga()

