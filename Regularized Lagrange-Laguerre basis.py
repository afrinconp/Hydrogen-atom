import matplotlib.pyplot as plt
from scipy.special import eval_laguerre
from scipy.special import roots_laguerre 
import numpy as np
from scipy import optimize

r=np.arange(0,20,0.001)

def root(N): #zeros laguerre
	a=roots_laguerre(N)[0]
	return a

def L(N,x): #eval laguerre function
	a=eval_laguerre (N,x)
	return a
	
def f(N,r):
	M=[]
	zerosN=root(N)
	k=0
	for xj in zerosN:
		p=(-1)**(k+1)*xj**(1./2.)*(r/xj)**n*(L(N,r)/(r-xj))*np.exp(-r/2)
		k=k+1
		M.append(p)
	return M

N=4
n=1
l=f(N,r)
plt.ylabel("fi(r)")
plt.xlabel("r")
plt.plot(r,l[0],label="i=1")
plt.plot(r,l[1],label="i=2")
plt.plot(r,l[2],label="i=3")
plt.plot(r,l[3],label="i=4")
plt.legend()
plt.show()
