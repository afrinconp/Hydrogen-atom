import matplotlib.pyplot as plt
from scipy.special import eval_legendre 
from scipy.special import roots_legendre  
import numpy as np
from scipy import optimize

r=np.arange(0,6,0.001)
a=6
def root(n):
	zero=roots_legendre(n)[0]
	xi=(zero+1)/2
	return xi
	


def f(N,r):
	x=2*r/a-1
	M=[]
	zerosN=root(N)
	k=0
	for xi in zerosN:
		print (xi)
		p=(-1)**(N+k+1)*(r/(a*xi))*(a*xi*(1-xi))**(1/2)*eval_legendre(N,x)/(r-a*xi)
		k=k+1
		M.append(p)
	return M
N=20
l=f(N,r)
plt.plot(r,l[2],label="i=3")
plt.plot(r,l[6],label="i=7")
plt.plot(r,l[9],label="i=10")
plt.plot(r,l[12],label="i=13")
plt.plot(r,l[16],label="i=17")
plt.legend()
plt.show()


