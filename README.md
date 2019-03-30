# Hydrogen-atom

## Regulazed Laguerre functions
The Regulazed Laguerre functions are defined in the interval (0,inf) as:


![Regulazed Laguerre functions](src/RegularizedLagrangeLaguerre.png)


The picture below shows the Regulazed Lagrange-Laguerre functions for N=4 and n=1. Basis functions f<sub>i</sub>(r) are plotted for i=1,2,3 and 4.


![Regulazed Laguerre functions](RegularizedLagrange-LaguerreFunctions.png)

L<sub>N</sub>(r) is a Laguerre polynomial of order N and the mesh points x<sub>j</sub> are the zeros of this polynomial. The kinetic-matrix elements at the Gauss approximation read

![LaguerreKineticMatrix1](src/LaguerreKineticMatrix1.png)

![LaguerreKineticMatrix2](src/LaguerreKineticMatrix2.png)

The Schrödinger equation can be written as:

![SchroEquation](src/SchroEquation.png)


It was decided to take u=xh and expand the wave function on the scaled-Lagrange basis defined in the interval (ha,hb) as

![phiu](src/phiu.png)

Then the wave equation is inserted in the Schrödinger equation and taking the properties of Gauss Quadrature, we end up with :

![PropiosVectors](src/PropiosVectors.png)


The expression above is a matrix so we can get information from it if we diagonalizate it. We can get eigenvalues which are the energies and the respecting eigenvector. After that, we can construct the wave function.

I constructed a program in python that diagonalizate the matrix, find thw wave function and then plot the probability density function (Laguerre-probability-density-function.py).



