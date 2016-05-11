import numpy as np

def parameters(data):
	T = 8
	dt = 0.001
	file = open(data, 'r')
	cd = []
	cl = []
	delP = []
	for i in range(int(T/dt)-20):
		line = file.readline()
		word = line.split()
		cd.append(float(word[0]))
		cl.append(float(word[1]))
		delP.append(float(word[2]))
		
 	maxcl = max(cl)
 	maxcd = max(cd)
 	t0 = []
 	for j in range(len(cl)):
 		if abs(cl[j] - maxcl) < 1e-4:
 			t0.append(j)

 	Period1 = 0
 	Period = np.zeros(2)
 	for k in range(len(t0)-1):
 		diff = abs(t0[k]-t0[k+1])/2
 		if diff > Period1:
 			Period[0] = t0[k]
 			Period[1] = t0[k+1]

 	t = [a for a in range(int(Period[0]),int(Period[1])+1)]
 	Tp = (t[-1]-t[0])
 	#f = lambda t: cl[t]

 	print maxcl, maxcd, 'St:', 0.1*10**3*(1./Tp)/(2*1.5/3.) # D*f/Ua
 	print 'delp:', delP[t[0]+int(round(0.5*Tp))]#np.sum(delP[t[0]:t[-1]+1])/len(delP[t[0]:t[-1]+1])

 	import matplotlib.pyplot as plt
 	plt.plot(t, cl[t[0]:t[-1]+1], label='Cl')
 	plt.plot(t, cd[t[0]:t[-1]+1], label='Cd')
 	plt.plot(t, delP[t[0]:t[-1]+1], label='deltaP')
 	plt.xlabel('t [s*10^-3]')
 	plt.legend(loc=4)
 	plt.title('One Period, starting from t at Clmax')
 	plt.show()



parameters('output1.txt')