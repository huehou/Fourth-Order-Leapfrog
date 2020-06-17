import Leapfrog4 as lf
import matplotlib.pyplot as plt
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
    # return np.abs(res - 1)

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

    # return overlap_2d(psi0, psi1, space), end-start

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

plt.show()
