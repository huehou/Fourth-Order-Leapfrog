import Leapfrog4 as lf
import matplotlib.pyplot as plt
import scipy.special as special
import numpy as np
import time

# LaTeX for plots
plt.rc('text', usetex = True)
plt.rc('font', family = 'serif', size = 18)

# Color scheme for plots
color = {
        'U72': '#00dcdc',
        'U3': '#0000ff',
        'U11': '#00ff00',
        'U7': '#ff0000',
        'RK4': '#ffa500'
        }

# ==============
# Pendulum
# ==============

# Parameters
m = 1 # mass
g = 9.81 # gravitational acceleration
omega = 1 # angular frequency
a = g/omega**2 # length of pendulum

# g = 9.81
# l = 1
# omega = np.sqrt(g/l)

# Potential energy
def V(x):
    return m * a**2 * omega**2 * (1 - np.cos(x))

def dV(x):
    return m * omega**2 * np.sin(x)

def dV2(x):
    return m * omega**2 * np.cos(x)

# Initial condition
t0 = 0
x0 = np.pi/3
p0 = 0

# Period
energy = V(x0)
period = 4/omega * special.ellipk(energy/(2*m*(omega**2)*(a**2)))

# Error
def error(xlist, plist):
    x1 = xlist[-1]
    p1 = plist[-1]
    return np.sqrt((x1 - x0)**2 + (p1 - p0)**2)

# Period Performance
def PeriodPerformance(option, x0, p0, dt, step, dV, dV2 = None):
    if option == 'U3':
        evolve = lf.U3C_1step
    elif option == 'U72':
        evolve = lf.U72C_1step
    elif option == 'RK4':
        evolve = lf.RK4C_1step
    elif option == 'U11':
        evolve = lf.U11C_1step
    elif option == 'U7':
        evolve = lf.U7C_1step
    else:
        raise ValueError("Option is not supported")
    
    x1, p1 = x0, p0

    start = time.time()
    if option != 'U7':
        for i in range(step):
            x1, p1 = evolve(x1, p1, dt, dV)
    else:
        for i in range(step):
            x1, p1 = evolve(x1, p1, dt, dV, dV2)


    end = time.time()

    return np.sqrt((x1 - x0)**2 + (p1 - p0)**2), end-start

steplist = np.logspace(1, 3, 100)
errorU3 = np.zeros(len(steplist))
errorU72 = np.zeros(len(steplist))
errorRK4 = np.zeros(len(steplist))
errorU11 = np.zeros(len(steplist))
errorU7 = np.zeros(len(steplist))
timeU3 = np.zeros(len(steplist))
timeU72 = np.zeros(len(steplist))
timeRK4 = np.zeros(len(steplist))
timeU11 = np.zeros(len(steplist))
timeU7 = np.zeros(len(steplist))

for i, step in enumerate(steplist):
    dt = period/int(step)
    errorU3[i], timeU3[i] = PeriodPerformance('U3', x0, p0, dt, int(step), dV)
    errorU72[i], timeU72[i] = PeriodPerformance('U72', x0, p0, dt, int(step), dV)
    errorRK4[i], timeRK4[i] = PeriodPerformance('RK4', x0, p0, dt, int(step), dV)
    errorU11[i], timeU11[i] = PeriodPerformance('U11', x0, p0, dt, int(step), dV)
    errorU7[i], timeU7[i] = PeriodPerformance('U7', x0, p0, dt, int(step), dV, dV2)

# Number of steps vs period error
plt.loglog(errorU3, steplist, '+', markersize = 6, label = "$U_3$", color = color["U3"])
plt.loglog(errorU72, steplist, '+', markersize = 6, label = "$U_7'$", color = color["U72"])
plt.loglog(errorRK4, steplist, '+', markersize = 6, label = "RK4", color = color["RK4"])
plt.loglog(errorU11, steplist, '+', markersize = 6, label = "$U_{11}$", color = color["U11"])
plt.loglog(errorU7, steplist, '+', markersize = 6, label = "$U_7$", color = color["U7"])
plt.legend(frameon = False)
plt.xlabel("period error $\epsilon_\mathrm{C}$")
plt.ylabel("number of steps N")
plt.tick_params(right = True, top = True, direction = 'in')

# Computation time vs period error
plt.figure()
plt.loglog(errorU3, timeU3, '+', markersize = 6, label = "$U_3$", color = color["U3"])
plt.loglog(errorU72, timeU72, '+', markersize = 6, label = "$U_7'$", color = color["U72"])
plt.loglog(errorRK4, timeRK4, '+', markersize = 6, label = "RK4", color = color["RK4"])
plt.loglog(errorU11, timeU11, '+', markersize = 6, label = "$U_{11}$", color = color["U11"])
plt.loglog(errorU7, timeU7, '+', markersize = 6, label = "$U_7$", color = color["U7"])
plt.legend(frameon = False)
plt.xlabel("period error $\epsilon_\mathrm{C}$")
plt.ylabel("computation time [s]")
plt.tick_params(right = True, top = True, direction = 'in')
plt.ylim(5*10**(-6), 7*10**(-1))


plt.show()
