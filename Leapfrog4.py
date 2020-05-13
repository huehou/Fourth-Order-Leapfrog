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

if __name__=="__main__":
    def V(x):
        return 1/2*np.array(x)**2

    def dV(x):
        return np.array(x)

    def dV2(x):
        return [[1, 0],[0, 1]]

    print(U7C_1step(dV, dV2, [1,1], [1,0], 0.1))
    print(RK4C_1step(dV, [1,1], [1,0], 0.1))
    print(U3C_1step(dV, [1,1], [1,0], 0.1))
    print(U72C_1step(dV, [1,1], [1,0], 0.1))
    print(U11C_1step(dV, [1,1], [1,0], 0.1))
    # print(RK4C_2d_1step(dV2, [1,1], [1,0], 0.1))
