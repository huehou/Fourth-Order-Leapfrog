import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import time
import scipy.interpolate as inter
from mpl_toolkits.mplot3d import Axes3D
import os

# =====================================
# Basic Functions - Classical
# =====================================

def U3C_1step(dV, x0, p0, dt):
    '''
    One step of standard leapfrog method for classical dynamics
    input: - dV: Function for first derivative of potential energy
           - x0: Starting position
           - p0: Starting momentum
           - dt: Size of timestep
    output: - x1: New position after one step
            - p1: New momentum after one step
    '''
    x0 = np.array(x0)
    p0 = np.array(p0)

    p_half = p0 - 0.5*dt*dV(x0)
    x1 = x0 + dt*p_half
    p1 = p_half - 0.5*dt*dV(x1)

    return x1, p1

def RK4C_1step(dV, x0, p0, dt):
    '''
    One step of 4th-order Runge-Kutta method for classical dynamics
    input: - dV: Function for first derivative of potential energy
           - x0: Starting position
           - p0: Starting momentum
           - dt: Size of timestep
    output: - x1: New position after one step
            - p1: New momentum after one step
    '''
    x0 = np.array(x0)
    p0 = np.array(p0)

    k1x = p0
    k1p = -dV(x0)
    
    k2x = p0 + dt/2*k1p
    k2p = - dV(x0 + dt/2*k1x)
    
    k3x = p0 + dt/2*k2p
    k3p = - dV(x0 + dt/2*k2x)
    
    k4x = p0 + dt*k3p
    k4p = - dV(x0 + dt*k3x)
    
    p1 = p0 + dt/6*(k1p + 2*k2p + 2*k3p + k4p)
    x1 = x0 + dt/6*(k1x + 2*k2x + 2*k3x + k4x)
    
    return x1, p1

def U7C_1step(dV, dV2, x0, p0, dt):
    '''
    One step of 4th-order leapfrog method for classical dynamics
    input: - dV: Function for first derivative of potential energy
           - dV2: Function for second derivative of potential energy
           - x0: Starting position
           - p0: Starting momentum
           - dt: Size of timestep
    output: - x1: New position after one step
            - p1: New momentum after one step
    '''
    x0 = np.array(x0)
    p0 = np.array(p0)

    p_one6 = p0 - dt/6*dV(x0)
    x_one2 = x0 + dt*p_one6/2
    p_five6 = p_one6 - dt*2/3*dV(x_one2) + dt**3/36*np.dot(dV2(x_one2),dV(x_one2))
    x1 = x_one2 + 0.5*p_five6*dt
    p1 = p_five6 - dt/6*dV(x1)
    
    return x1, p1

def U72C_1step(dV, x0, p0, dt):
    '''
    One step of U7' method for classical dynamics
    input: - dV: Function for first derivative of potential energy
           - x0: Starting position
           - p0: Starting momentum
           - dt: Size of timestep
    output: - x1: New position after one step
            - p1: New momentum after one step
    '''
    x0 = np.array(x0)
    p0 = np.array(p0)
    s = 1/(2 - 2**(1/3))

    pfirst = p0 - s/2*dt*dV(x0)
    xsecond = x0 + s*dt*pfirst
    pthird = pfirst - (1-s)/2*dt*dV(xsecond)
    xfourth = xsecond + (1-2*s)*dt*pthird
    pfifth = pthird - (1-s)/2*dt*dV(xfourth)
    x1 = xfourth + s*dt*pfifth
    p1 = pfifth - s/2*dt*dV(x1)

    return x1, p1

def U11C_1step(dV, x0, p0, dt):
    '''
    One step of U11 method for classical dynamics 
    input: - dV: Function for first derivative of potential energy
           - x0: Starting position
           - p0: Starting momentum
           - dt: Size of timestep
    output: - x1: New position after one step
            - p1: New momentum after one step
    '''
    s = 1/(4 - 4**(1/3))
    
    pfirst = p0 - s/2*dt*dV(x0)
    xsecond = x0 + s*dt*pfirst
    pthird = pfirst - s*dt*dV(xsecond)
    xfourth = xsecond + s*dt*pthird
    pfifth = pthird - (1-3*s)/2*dt*dV(xfourth)
    xsixth = xfourth + (1-4*s)*dt*pfifth
    pseventh = pfifth - (1-3*s)/2*dt*dV(xsixth)
    xeighth = xsixth + s*dt*pseventh
    pninth = pseventh - s*dt*dV(xeighth)
    xtenth = xeighth + s*dt*pninth
    peleventh = pninth - s/2*dt*dV(xtenth)

    x1, p1 = xtenth, peleventh

    return x1, p1

def evolveC(option, t0, x0, p0, dt, steps, dV, dV2=None):
    '''
    Evolution for classical dynamics
    input: - option: U3, RK4, U7, U72, or U11 algorithm
           - t0: Starting time
           - x0: Starting position
           - p0: Starting momentum
           - dt: Size of timestep
           - steps: Number of steps forward
           - dV: Function for first derivative of potential energy
           - dV2: Function for second derivative of potential energy (optional)
    output: - tlist: List of time
            - xlist: List of position
            - plist: List of momentum
    '''
    if option=="U7" and dV2==None:
        raise ValueError("Second derivative is needed for U7")
    if type(x0) == type(1) or type(x0) == type(1.):
        dim = 1
    else:
        dim = len(x0) # Dimension of the system

    xlist = np.zeros((steps+1, dim))
    plist = np.zeros((steps+1, dim))
    tlist = np.zeros(steps+1)
    
    xlist[0] = np.array(x0)
    plist[0] = np.array(p0)
    tlist[0] = t0

    for i in range(1, steps+1):
        if option == "U3":
            x1, p1 = U3C_1step(dV, xlist[i-1], plist[i-1], dt)
        elif option == "RK4":
            x1, p1 = RK4C_1step(dV, xlist[i-1], plist[i-1], dt)
        elif option == "U7":
            x1, p1 = U7C_1step(dV, dV2, xlist[i-1], plist[i-1], dt)
        elif option == "U72":
            x1, p1 = U72C_1step(dV, xlist[i-1], plist[i-1], dt)
        elif option == "U11":
            x1, p1 = U11C_1step(dV, xlist[i-1], plist[i-1], dt)
        else:
            raise ValueError("Option is not valid")
        xlist[i] = x1
        plist[i] = p1
        tlist[i] = tlist[i-1] + dt

    return np.array(tlist), np.array(xlist), np.array(plist)

# ========================================
# Basic Functions - Fourier Transforms
# ========================================

def fft(x, fun)
    '''
    Perform fast Fourier transform with corrected factors in 1D
    input: - x: Position coordinates
           - fun: Function in position space
    output: - k: Wave number coordinates
            - fun: Function in wave number space
    '''
    
    # Normal FFT procedures
    fun = np.fft.ifftshift(fun)
    fun = np.fft.fft(fun)
    fun = np.fft.fftshift(fun)

    # Multiply by correcting factors
    fun = fun * (x[1] - x[0]) / np.sqrt(2*np.pi)

    # Wave number space
    k = np.fft.fftfreq(len(x), x[1]-x[0])
    k = np.fft.fftshift(k)
    k = k * 2 * np.pi

    return k, fun

def ifft(k, fun):
    '''
    Perform inverse fast Fourier transform with corrected factors in 1D
    input: - k: Wave number coordinates
           - fun: Function in wave number space
    output: - x: Position coordinates
            - fun: Function in position space
    '''

    # Normal IFFT procedures
    fun = np.fft.ifftshift(fun)
    fun = np.fft.ifft(fun)
    fun = np.fft.fftshift(fun)

    # Multiply by correcting factors
    fun = fun * (k[1] - k[0]) / np.sqrt(2 * np.pi) * len(k)

    # Wave number space
    x = np.fft.fftfreq(len(k), k[1]-k[0])
    x = np.fft.fftshift(x)
    x = x * 2 * np.pi

    return x, fun


if __name__=="__main__":
    def V(x):
        return 1/2*np.array(x)**2

    def dV(x):
        return np.array(x)

    def dV2(x):
        # return 1
        return [[1, 0],[0, 1]]

    x0 = [1,1]
    p0 = [-1,0]
    resU3 = evolveC("U3", 0, x0, p0, 0.01, 100, dV)
    resRK4 = evolveC("RK4", 0, x0, p0, 0.01, 100, dV)
    resU7 = evolveC("U7", 0, x0, p0, 0.01, 100, dV, dV2)
    plt.plot(resU3[1][:,0],resU3[1][:,1])
    plt.plot(resRK4[1][:,0],resRK4[1][:,1])
    plt.plot(resU7[1][:,0],resU7[1][:,1])
    plt.show()
