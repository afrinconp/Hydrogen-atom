#This program gives the laguerre wave function base state for hidrogen atom 
import matplotlib.pyplot as plt
from scipy.special import eval_laguerre 
from scipy.special import roots_laguerre 
import scipy.integrate as integrate
import numpy as np

#-------------------------------Laguerre roots and function--------------------------------#

 
def root(x): #zeros laguerre
	a=roots_laguerre(x)[0]
	return a

def L(N,x): #eval laguerre function
	a=eval_laguerre (N,x)
	return a

#-----------------------------------------Parameters---------------------------------------#
	
e=2      #Costante electron  
N=30
n=1
h=1.5

#----------------------------------------Matrix--------------------------------------------#

def v(x):#Coulumb potential
	a=-e/x
	return a
	
zeros1=root(N)
zeros2=root(N)

def T(zeros1,zeros2,n,N): #Matrix Elements
	a = np.matrix(np.zeros((N, N), dtype = np.float))
	z=0
	for xi in zeros2:
		k=0
		for xj in zeros1:
			if xi!=xj:
				element=(1/h**2)*(-1)**(-k+z+1)*(xj**(n-3./2.)/xi**(n-1./2))*((2*n-3)*xj-(2*n-1)*xi)/(xi-xj)**2
				a[z,k]=element
			elif xi==xj:
				a[z,k]=(1/h**2)*(12*xi**2)**(-1)*(-12*n**2+24*n-8+(4*N+2)*xi-xi**2)+v(xi*h)
			k=k+1
		z=z+1
	return a

c=T(zeros1,zeros2,n,N)
#----------------------------------Vectors and energies from the matrix-----------------------#

d,vectors=np.linalg.eig(c)                #Diagonaliza la matriz
eigenvalues=np.sort(d)                    #Sort Energies
Energybase=eigenvalues[0]                 #Energy base state
listd=np.argwhere(d==Energybase)          #It finds where the energy base state is in d
p=listd.tolist()                          #Conver numpy to alist 
position=p[0][0]                          #The position of energy base state in vectors 
vectorbaseA=vectors[:,int(position)]      #Vector of enegy base state in matrix form
vectorbaseAA=vectorbaseA.tolist()         #Convert numpy matrix vectorbase to a list
k=[]                                      #Create a list
for i in vectorbaseAA:
	k.append(i[0])
vectorbase=np.array(k)                    #Vector numpy array base state

print ("Vector - estado base",vectorbase)

#-------------------------------Laguerre Regulazed function-------------------------------------#

r=np.arange(0,25,0.001)
def f(N,r):
	M=[]
	zerosN=root(N)
	k=0
	for xj in zerosN:
		p=((-1)**(k+1)*xj**(1./2.)*(r/(h*xj))**n*(L(N,r/h)/((r/h)-xj))*np.exp(-r/(h*2.)))/(h**(1./2))
		k=k+1
		M.append(p)
	return M	 
Regulazed_Laguerre_functions=f(N,r)



#--------------------------------Wave function--------------------------------------------------#
def Wave(x):
	k=0
	Wa=[]
	for i in vectorbase:
		a=i*Regulazed_Laguerre_functions[k]  #Components of vector base multiplied by fuction
		k=k+1
		Wa.append(a)
	W=Wa[0]
	c=0
	for i in range(len(Wa)-1):               #Plus componentes of Wa
		c=c+1
		W=W+Wa[c]
	return W
	
phi=Wave(r)                       #Wave function
prob=phi*phi                      #Probability Density
		 

coeficientes2=np.polyfit(r,prob,15)        #Ajuste polinomico
polinomio2=np.poly1d(coeficientes2)


int2=integrate.quad(polinomio2, 0, 25)     #Integral
print ("integral probabilidad", int2)

plt.plot(r,prob,label="Probability density function")
plt.legend()
plt.show()


		

	
