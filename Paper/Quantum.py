import Leapfrog4 as lf
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as special
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


# # ===========================
# # 2D Harmonic Osciilator
# # ===========================

# # Potential energy
# def V(X, Y):
    # return 1/2*(X**2 + Y**2)

# def dV(X, Y):
    # return np.array([X, Y])

# # Initial conditions
# t0 = 0
# dt = 0.1
# x = np.arange(-10, 10, 1/20)
# y = np.arange(-10, 10, 1/20)
# X, Y = np.meshgrid(x, y, indexing = 'ij')
# space = np.array([X, Y])
# psi0 = (X - 1j*Y)*np.exp(-((X-1)**2)/2) * np.exp(-((Y-1)**2)/2) / np.sqrt(3*np.pi)
# period = 2*np.pi

# # Overlap between two wave function
# def overlap_2d(psi1, psi2, space):
    # product = np.conj(psi1)*psi2
    # res = np.trapz(np.trapz(product, space[0][:,0]), space[1][0,:])
    # return res

# # Period Performance
# def PeriodPerformance(option, space, psi0, dt, step, V, dV = None):
    # if option == 'U3':
        # evolve = lf.U3Q_2d_1step
    # elif option == 'U72':
        # evolve = lf.U72Q_2d_1step
    # elif option == 'RK4':
        # evolve = lf.RK4Q_2d_1step
    # elif option == 'U11':
        # evolve = lf.U11Q_2d_1step
    # elif option == 'U7':
        # evolve = lf.U7Q_2d_1step
    # else:
        # raise ValueError("Option is not supported")
    
    # psi1 = psi0

    # start = time.time()
    # if option != 'U7':
        # for i in range(step):
            # space, psi1 = evolve(space, psi1, dt, V)
    # else:
        # for i in range(step):
            # space, psi1 = evolve(space, psi1, dt, V, dV)


    # end = time.time()

    # return np.abs(overlap_2d(psi0, psi1, space)-1), end-start

# steplist = np.logspace(1, 2.5, 100)
# errorU3 = np.zeros(len(steplist))
# errorU72 = np.zeros(len(steplist))
# errorU11 = np.zeros(len(steplist))
# errorU7 = np.zeros(len(steplist))
# timeU3 = np.zeros(len(steplist))
# timeU72 = np.zeros(len(steplist))
# timeU11 = np.zeros(len(steplist))
# timeU7 = np.zeros(len(steplist))

# for i, step in enumerate(steplist):
    # dt = period/int(step)
    # print("Progress:", int(i/len(steplist)*100))
    # errorU3[i], timeU3[i] = PeriodPerformance("U3", space, psi0, dt, int(step), V)
    # errorU72[i], timeU72[i] = PeriodPerformance("U72", space, psi0, dt, int(step), V)
    # errorU11[i], timeU11[i] = PeriodPerformance("U11", space, psi0, dt, int(step), V)
    # errorU7[i], timeU7[i] = PeriodPerformance("U7", space, psi0, dt, int(step), V, dV)

# plt.loglog(errorU3, steplist, '.', markersize = 2, color = color["U3"], label = "$U_3$")
# plt.loglog(errorU72, steplist, '.', markersize = 2, color = color["U72"], label = "$U_7'$")
# plt.loglog(errorU11, steplist, '.', markersize = 2, color = color["U11"], label = "$U_{11}$")
# plt.loglog(errorU7, steplist, '.', markersize = 2, color = color["U7"], label = "$U_7$")
# plt.tick_params(direction = "in", top = True, right = True)
# plt.legend(frameon = False)
# plt.xlabel("period error $\epsilon_\mathrm{Q}$")
# plt.ylabel("number of steps $N$")

# plt.figure()
# plt.loglog(errorU3, timeU3, '.', markersize = 2, color = color["U3"], label = "$U_3$")
# plt.loglog(errorU72, timeU72, '.', markersize = 2, color = color["U72"], label = "$U_7'$")
# plt.loglog(errorU11, timeU11, '.', markersize = 2, color = color["U11"], label = "$U_{11}$")
# plt.loglog(errorU7, timeU7, '.', markersize = 2, color = color["U7"], label = "$U_7$")
# plt.tick_params(direction = "in", top = True, right = True)
# plt.legend(frameon = False)
# plt.xlabel("period error $\epsilon_\mathrm{Q}$")
# plt.ylabel("computation time [s]")

# =====================
# Davidson Potential
# =====================

# Parameters
m = 1
omega = 1
hbar = 1

# Potential energy
def V(X, Y, Z):
    res = 1/2*m*omega**2 * (X**2 + Y**2 + Z**2) + np.true_divide(hbar**2, (2*m*(X**2 + Y**2 + Z**2)))
    res[res == np.inf] = 0
    res = np.nan_to_num(res)
    return res

def dV(X, Y, Z):
    dVx = m*omega**2 * X - (X*hbar)/(m*(X**2 + Y**2 + Z**2)**2)
    dVy = m*omega**2 * Y - (Y*hbar)/(m*(X**2 + Y**2 + Z**2)**2)
    dVz = m*omega**2 * Z - (Z*hbar)/(m*(X**2 + Y**2 + Z**2)**2)
    
    return np.array([dVx, dVy, dVz])

x = np.linspace(-10, 10, 101)
y = x
z = x
X, Y, Z = np.meshgrid(x, y, z, indexing = 'ij')
space = [X, Y, Z]

# Davidson potential eigenstates
def basis(n, l, m, X, Y, Z):
    R = np.sqrt(X**2 + Y**2 + Z**2)
    Theta = np.arctan2(np.sqrt(X**2 + Y**2), Z)
    Phi = np.arctan2(Y, X)
    alpha = 1/4*(np.sqrt(4*l*(l+1) + 4 + 1) - 1)
    beta = 1/2
    gamma = 2*alpha + 1/2

    res = np.exp(-R**2/2) * (R**2)**alpha + 0*1j
    res *= special.sph_harm(m, l, Phi, Theta)
    res *= special.eval_genlaguerre(n, 2*alpha + 1/2, R**2)

    return res

# Exact time evolution of Davidson bases
def exact(n, l, m, X, Y, Z, t):
    res = basis(n, l, m, X, Y, Z)
    alpha = 1/4*(np.sqrt(4*l*(l+1) + 4 + 1) - 1)
    energy = 2*n + 3/2 + 2*alpha
    res *= np.exp(-1j*t*energy)

    return res

# Overlap between two wave function
def overlap_3d(psi1, psi2, space):
    product = np.conj(psi1)*psi2
    # res = np.trapz(np.trapz(np.trapz(np.conj(psi1)*psi2, space[0][:,0][:,0]), space[1][0,:][:,0]), space[2][0,:][0,:])
    res = np.trapz(np.trapz(np.trapz(product, space[0][:,0][:,0]), space[1][0,:][:,0]), space[2][0,:][0,:])
    return res

# Initial state
n = 1
l = 20
m = 20
psi0 = basis(n, l, m, X, Y, Z)
norm = np.abs(overlap_3d(psi0, psi0, space)) # normalise wave function
psi0 /= np.sqrt(norm)

alpha = 1/4*(np.sqrt(4*l*(l+1) + 4 + 1) - 1)
energy = 2*n + 3/2 + 2*alpha
period = 2*np.pi/energy

# Period performance
def PeriodPerformance(option, space, psi0, dt, step, V, dV = None):
    if option == 'U3':
        evolve = lf.U3Q_3d_1step
    elif option == 'U72':
        evolve = lf.U72Q_3d_1step
    elif option == 'RK4':
        evolve = lf.RK4Q_3d_1step
    elif option == 'U11':
        evolve = lf.U11Q_3d_1step
    elif option == 'U7':
        evolve = lf.U7Q_3d_1step
    else:
        raise ValueError("Option is not supported")
    
    psi1 = psi0

    start = time.time()
    if option != 'U7':
        for i in range(step):
            space, psi1 = evolve(space, psi1, dt, V)
    else:
        for i in range(step):
            space, psi1 = evolve(space, psi1, dt, V, dV)


    end = time.time()

    return np.abs(overlap_3d(psi0, psi1, space)-1), end-start

steplist = np.logspace(0, 2, 100)
errorU3 = np.zeros(len(steplist))
errorU72 = np.zeros(len(steplist))
errorU11 = np.zeros(len(steplist))
errorU7 = np.zeros(len(steplist))
timeU3 = np.zeros(len(steplist))
timeU72 = np.zeros(len(steplist))
timeU11 = np.zeros(len(steplist))
timeU7 = np.zeros(len(steplist))

for i, step in enumerate(steplist):
    dt = period/int(step)
    print("Progress:", int(i/len(steplist)*100))
    errorU3[i], timeU3[i] = PeriodPerformance("U3", space, psi0, dt, int(step), V)
    errorU72[i], timeU72[i] = PeriodPerformance("U72", space, psi0, dt, int(step), V)
    errorU11[i], timeU11[i] = PeriodPerformance("U11", space, psi0, dt, int(step), V)
    errorU7[i], timeU7[i] = PeriodPerformance("U7", space, psi0, dt, int(step), V, dV)

plt.loglog(errorU3, steplist, '.', markersize = 2, color = color["U3"], label = "$U_3$")
plt.loglog(errorU72, steplist, '.', markersize = 2, color = color["U72"], label = "$U_7'$")
plt.loglog(errorU11, steplist, '.', markersize = 2, color = color["U11"], label = "$U_{11}$")
plt.loglog(errorU7, steplist, '.', markersize = 2, color = color["U7"], label = "$U_7$")
plt.tick_params(direction = "in", top = True, right = True)
plt.legend(frameon = False)
plt.xlabel("period error $\epsilon_\mathrm{Q}$")
plt.ylabel("number of steps $N$")

plt.figure()
plt.loglog(errorU3, timeU3, '.', markersize = 2, color = color["U3"], label = "$U_3$")
plt.loglog(errorU72, timeU72, '.', markersize = 2, color = color["U72"], label = "$U_7'$")
plt.loglog(errorU11, timeU11, '.', markersize = 2, color = color["U11"], label = "$U_{11}$")
plt.loglog(errorU7, timeU7, '.', markersize = 2, color = color["U7"], label = "$U_7$")
plt.tick_params(direction = "in", top = True, right = True)
plt.legend(frameon = False)
plt.xlabel("period error $\epsilon_\mathrm{Q}$")
plt.ylabel("computation time [s]")



plt.show()
