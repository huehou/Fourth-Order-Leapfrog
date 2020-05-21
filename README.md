# Fourth-Order-Leapfrog

This is an online open-source code complimentary to the article [Fourth-order leapfrog algorithms for numerical time evolution of classical and quantum systems](). All algorithms in the paper are defined in Leapfrog4.py. To use the code, simply copy Leapfrog4.py to the repository and import it into python 

~~~~
import Leapfrog4 as lf
~~~~

Here is a summary of useful functions in the module.

## Classical Dynamics

The functions for classical evolutions are
- U3C: classical algorithm for $U_3$ or the common leapfrog algorithm
- RK4C: classical algorithm for 4th-order Runge--Kutta
- U72C: classical algorithm for $U_7'$, one of the 4th-order leapfrog algorithms
- U11C: classical algorithm for $U_{11}$, one of the 4th-order leapfrog algorithms
- U7C: classical algorithm for $U_7$, one of the 4th-order leapfrog algorithms

They take inputs
- t0: the initial time
- x0: the initial position
- p0: the initial momentum
- dt: the size of time step
- steps: the number of steps
- dV: the function for the first derivative of the potential energy
- dV2: (only for U7C) the function for the second derivative of the potential energy  <br>

and give outputs
- tlist: a list of time
- xlist: a list of position
- plist: a list of momentum

### 1D Example

For example, consider classical harmonic oscillator in 1D. First define its first derivative dV, and the second derivative of the potential energy dV2.

~~~~
import numpy as np
import matplotlib.pyplot as plt

def dV(x):
	return x

def dV2(x):
	return 1
~~~~

Next, define the initial time t0, initial position x0, and initial momentum

~~~~
t0 = 0
x0 = 1
p0 = 0
~~~~

Then, decide the size of time step and the number of steps needed for the evolution. As an example, since the period is known for the 1D harmonic oscillator potential

~~~~
period = 2*np.pi
steps = 1000
dt = period/steps*3
~~~~

Finally, pick one evolution method and start the simulation. For illustration purposes, we will run all 5 functions here

~~~~
t1U3, x1U3, p1U3 = lf.U3C(t0, x0, p0, dt, steps, dV)
t1RK4, x1RK4, p1RK4 = lf.RK4C(t0, x0, p0, dt, steps, dV)
t1U72, x1U72, p1U72 = lf.U72C(t0, x0, p0, dt, steps, dV)
t1U11, x1U11, p1U11 = lf.U11C(t0, x0, p0, dt, steps, dV)
t1U7, x1U7, p1U7 = lf.U7C(t0, x0, p0, dt, steps, dV, dV2)
~~~~

To plot the result from $U_3$, for example, simply

~~~~
plt.plot(t1U3, x1U3)
~~~~

### 2D Example

Similarly, consider the example of 2D harmonic oscillator. The derivatives of the potential energies are 

~~~~
def dV(x):
	return np.array(x)

def dV2(x):
	return np.array([[1,0], [0,1]])
~~~~

For numerical convenience, the `numpy` array should be used. Furthermore, in 2D, dV is the gradient vector, while dV2 is the Hessian matrix.

The rest of the procedure follows similarly as before

~~~~
t0 = 0
x0 = np.array([0, 1])
p0 = np.array([10, 0])
period = 2*np.pi
steps = 100
dt = period/steps

t1U3, x1U3, p1U3 = lf.U3C(t0, x0, p0, dt, steps, dV)
t1RK4, x1RK4, p1RK4 = lf.RK4C(t0, x0, p0, dt, steps, dV)
t1U72, x1U72, p1U72 = lf.U72C(t0, x0, p0, dt, steps, dV)
t1U11, x1U11, p1U11 = lf.U11C(t0, x0, p0, dt, steps, dV)
t1U7, x1U7, p1U7 = lf.U7C(t0, x0, p0, dt, steps, dV, dV2)
~~~~

### 3D Example

The exact same procedure applies to 3D evolution

~~~~
def dV(x):
    return np.array(x)

def dV2(x):
    return np.array([[1,0,0],[0,1,0],[0,0,1]])

t0 = 0
x0 =np.array([0,1,0])
p0 = np.array([1,0,1])
period = 2*np.pi
steps = 100
dt = period/steps

t1U3, x1U3, p1U3 = lf.U3C(t0, x0, p0, dt, steps, dV)
t1RK4, x1RK4, p1RK4 = lf.RK4C(t0, x0, p0, dt, steps, dV)
t1U72, x1U72, p1U72 = lf.U72C(t0, x0, p0, dt, steps, dV)
t1U11, x1U11, p1U11 = lf.U11C(t0, x0, p0, dt, steps, dV)
t1U7, x1U7, p1U7 = lf.U7C(t0, x0, p0, dt, steps, dV, dV2)
~~~~

The examples for classical dynamics can be found in EgClassical.py.

### Misc

Leapfrog4.py also provides the algorithms 1 step at a time. Simply evolve them by
~~~~
x1, p1 = U3C_1step(x0, p0, dt, dV)
x1, p1 = RK4C_1step(x0, p0, dt, dV)
x1, p1 = U72C_1step(x0, p0, dt, dV)
x1, p1 = U11C_1step(x0, p0, dt, dV)
x1, p1 = U7C_1step(x0, p0, dt, dV, dV2)
~~~~
where x1 and p1 are the position and momentum after one time step, respectively. To generalise the program to include different algorithms, the program also defined a general evolution function `evolveC(option, t0, x0, p0, dt, steps, dV, dV2=None)`. Simply define a function for the new algorithm evolving for 1 time step, and add a respective option in the `evolveC` program.

## Quantum Dynamics

The functions for quantum evolutions are
- U3Q: quantum algorithm for $U_3$, quantum version of the common leapfrog algorithm
- U72Q: quantum algorithm for $U_7'$, one of the 4th-order unitary algorithms
- U11Q: quantum algorithm for $U_{11}$, one of the 4th-order unitary algorithms
- U7Q: quantum algorithm for $U_7$, one of the 4th-order unitary algorithms

They take inputs
- t0: initial time
- x0: position space
- psi0: initial wave function
- dt: the size of time step
- steps: the number of steps
- V: the function for potential energy
- dV: (only for U7Q) the first derivative of the potential energy

and give outputs

- tlist: a list of time
- x: the position space
- psilist: a list of wave functions

Since the algorithms involve Fourier transforms of the wave function. functions for (inverse) fast Fourier transforms with correction factors are provided

- fft: 1D fast Fourier transform
- ifft: 1D inverse fast Fourier transform
- fft: 2D fast Fourier transform
- ifft: 2D inverse fast Fourier transform
- fft: 2D fast Fourier transform
- ifft: 2D inverse fast Fourier transform

The (inverse) Fourier transform functions take inputs

- X (K): meshgrid of position (wave number) coordinates
- psi (phi): wave function in position (wave number) space

and give outputs

- K (X): meshgrid of wave number (position) coordinates
- phi (psi): wave function in wave number (position) space

### 1D Example

Consider 1D quantum harmonic oscillator potential

~~~~
def V(x):
	return 1/2*x**2

def dV(x):
	return x
~~~~

Next, define the position coordinates x and the initial wave function psi0.

~~~~
x = np.linspace(-10, 10, 1000)
psi0 = np.exp(-((x-1)**2/2))/(np.pi)**(1/4)
~~~~

Define the number of steps and size of time step
~~~~
t0 = 0
period = 4*np.pi
steps = 100
dt = period/steps
~~~~

Then simply apply the evolution functions to obtain a list of result
~~~~
t, x, psiU3 = lf.U3Q(t0, x, psi0, dt, steps, V)
t, x, psiU72 = lf.U72Q(t0, x, psi0, dt, steps, V)
t, x, psiU11 = lf.U11Q(t0, x, psi0, dt, steps, V)
t, x, psiU7 = lf.U7Q(t0, x, psi0, dt, steps, V, dV)
~~~~

To plot the final probability density of, say, the result from $U_7$ algorithm, simply plot
~~~~
plt.plot(x, np.abs(psiU7[-1])**2)
~~~~

If we are interested in the wave number space probability density, then the wave number wave function can be obtained by
~~~~
k, phiU7 = lf.fft(x, psiU7[-1])
~~~~

### 2D Example

Consider the 2D harmonic oscillator, As with the classical case, the first derivative of the potential energy is the gradient vector

~~~~
def V(X, Y):
	return 1/2*(X**2 + Y**2)

def dV(X, Y):
	return np.array([X, Y])
~~~~

We then define the meshgrid for the position space

~~~~
x = np.linspace(-10, 10, 100)
y = np.linspace(-10, 10, 100)
X, Y = np.meshgrid(x, y, indexing = "ij")
~~~~

Then the rest follows from before
~~~~
psi0 = np.exp(-(((X-1)**2+(Y-1)**2)/2))/(np.pi)**(1/2)
t0 = 0
period = 2*np.pi
steps = 50
dt = period/steps
~~~~

Then the evolution is done by
~~~~
t, S, psiU3 = lf.U3Q(t0, [X,Y], psi0, dt, steps, V)
t, S, psiU72 = lf.U72Q(t0, [X,Y], psi0, dt, steps, V)
t, S, psiU11 = lf.U11Q(t0, [X,Y], psi0, dt, steps, V)
t, S, psiU7 = lf.U7Q(t0, [X,Y], psi0, dt, steps, V, dV)
~~~~

The procedure for Fourier transform is 
~~~~
[kX, kY], phi = lf.fft2([X,Y], psi)
~~~~

### 3D Example

The same procedures follow from before
~~~~
def V(X, Y, Z):
    return 1/2*(X**2+Y**2+Z**2)

def dV(X, Y, Z):
    return np.array([X, Y, Z])

x = np.linspace(-3, 3, 50)
y = np.linspace(-3, 3, 50)
z = np.linspace(-3, 3, 50)
X, Y, Z = np.meshgrid(x, y, z, indexing = "ij")
psi0 = np.exp(-(((X-1)**2+(Y-1)**2+(Z-1)**2))/2)/(np.pi)**(3/4)
t0 = 0
period = 4*np.pi
steps = 50
dt = period/steps

t, S, psiU3 = lf.U3Q(t0, [X,Y,Z], psi0, dt, steps, V)
t, S, psiU72 = lf.U72Q(t0, [X,Y,Z], psi0, dt, steps, V)
t, S, psiU11 = lf.U11Q(t0, [X,Y,Z], psi0, dt, steps, V)
t, S, psiU7 = lf.U7Q(t0, [X,Y,Z], psi0, dt, steps, V, dV)
~~~~

Since 3D wave functions can get unreasonably large, it is sometimes impractical to store the entire list of wave function in the temporary memory. The program provides a method to save the wave functions into a temporary folder. This, using U7Q as an example, can be invoked by 

~~~~
t, S, psiU7 = lf.U7Q(t0, [X,Y,Z], psi0, dt, steps, V, dV, True)
~~~~

Once running, the program will prompt a choice of temporary directory to store the wave functions

~~~~
Input directory to store wave functions (default is ./temp):
~~~~

If you press `ENTER` without any arguments, the wave functions will either create a directory or save in the already existing directory "./temp". The output `psiU7` in this case is then a list of directories pointing to the respective wave functions as the corresponding time steps in `t`.

The examples discussed in this section can be found in EgQuantum.py.

### Misc

Similarly to the classical equivalent, Leapfrog4.py also provides the algorithms 1 step at a time. However, the 1-step evolution functions are dimension specific. Therefore, one can evolve U7, for example, 1 step at a time by 
~~~~
X, psi1 = lf.U7Q_1d_1step(X, psi0, dt, V, dV) # For 1D
[X,Y], psi1 = lf.U7Q_2d_1step([X,Y], psi0, dt, V, dV) # For 2D
[X,Y,Z], psi1 = lf.U7Q_3d_1step([X,Y,Z], psi0, dt, V, dV) # For 3D
~~~~

To generalise the program to include different algorithms, the program also defined a general evolution function `evolveQ_1d(option, t0, X, psi0, dt, steps, V, dV=None)`, `evolveQ_2d(option, t0, X, psi0, dt, steps, V, dV=None)`, and `evolveQ_3d(option, t0, X, psi0, dt, steps, V, dV=None)`. Simply define a function for the new algorithm evolving for 1 time step for the respective spatial dimensions, and add a respective option in the `evolveQ` programs
