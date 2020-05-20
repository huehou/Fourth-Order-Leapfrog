import numpy as np
import Leapfrog4 as lf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation


# # ===============================
# # Quantum 1D Harmonic Oscillator
# # ===============================

# def V(x):
    # return 1/2*x**2

# def dV(x):
    # return x

# def animate(lst, ax):
    # def func(i):
        # ax.clear()
        # line = ax.plot(x, np.abs(lst[i])**2)
        # return line,
    # return func

# x = np.linspace(-10, 10, 1000)
# psi0 = np.exp(-((x-1)**2/2))/(np.pi)**(1/4)
# t0 = 0
# period = 4*np.pi
# steps = 100
# dt = period/steps

# t, x, psiU3 = lf.U3Q(t0, x, psi0, dt, steps, V)
# t, x, psiU72 = lf.U72Q(t0, x, psi0, dt, steps, V)
# t, x, psiU11 = lf.U11Q(t0, x, psi0, dt, steps, V)
# t, x, psiU7 = lf.U7Q(t0, x, psi0, dt, steps, V, dV)

# def overlap(psi1, psi2, space):
    # res = np.trapz(np.conj(psi1)*psi2, space)
    # return res

# print("Error for U3 is {}".format(np.abs(overlap(psi0, psiU3[-1], x)-1)))
# print("Error for U72 is {}".format(np.abs(overlap(psi0, psiU72[-1], x)-1)))
# print("Error for U11 is {}".format(np.abs(overlap(psi0, psiU11[-1], x)-1)))
# print("Error for U7 is {}".format(np.abs(overlap(psi0, psiU7[-1], x)-1)))

# pltlst = psiU7
# fig = plt.figure()
# ax = fig.gca()
# ax.plot(x, np.abs(pltlst[0])**2)
# ani = animation.FuncAnimation(fig, animate(pltlst, ax), frames=steps, interval = 1, blit = False)

# # ===============================
# # Quantum 2D Harmonic Oscillator
# # ===============================

# def V(X, Y):
    # return 1/2*(X**2+Y**2)

# def dV(X, Y):
    # return np.array([X, Y])

# def animate(lst, ax):
    # def func(i):
        # ax.clear()
        # line = ax.contourf(S[0], S[1], np.abs(lst[i])**2)
        # return line,
    # return func

# x = np.linspace(-10, 10, 100)
# y = np.linspace(-10, 10, 100)
# X, Y = np.meshgrid(x, y, indexing = "ij")
# psi0 = np.exp(-(((X-1)**2+(Y-1)**2)/2))/(np.pi)**(1/2)
# t0 = 0
# period = 2*np.pi
# steps = 50
# dt = period/steps

# t, S, psiU3 = lf.U3Q(t0, [X,Y], psi0, dt, steps, V)
# t, S, psiU72 = lf.U72Q(t0, [X,Y], psi0, dt, steps, V)
# t, S, psiU11 = lf.U11Q(t0, [X,Y], psi0, dt, steps, V)
# t, S, psiU7 = lf.U7Q(t0, [X,Y], psi0, dt, steps, V, dV)

# def overlap(psi1, psi2, space):
    # res = np.trapz(np.trapz(np.conj(psi1)*psi2, space[0][:,0]), space[1][0,:])
    # return res

# print("Error for U3 is {}".format(np.abs(overlap(psi0, psiU3[-1], S)-1)))
# print("Error for U72 is {}".format(np.abs(overlap(psi0, psiU72[-1], S)-1)))
# print("Error for U11 is {}".format(np.abs(overlap(psi0, psiU11[-1], S)-1)))
# print("Error for U7 is {}".format(np.abs(overlap(psi0, psiU7[-1], S)-1)))

# pltlst = psiU7
# fig = plt.figure()
# ax = fig.gca()
# ax.contourf(S[0], S[1], np.abs(pltlst[0])**2)
# ani = animation.FuncAnimation(fig, animate(pltlst, ax), frames=steps, interval = 1, blit = False)

# ===============================
# Quantum 3D Harmonic Oscillator
# ===============================

def V(X, Y, Z):
    return 1/2*(X**2+Y**2+Z**2)

def dV(X, Y, Z):
    return np.array([X, Y, Z])

def animate(lst, ax):
    def func(i):
        ax.clear()
        ind = np.where(np.abs(pltlst[i])**2 > 0.05)
        lst = (S[0][ind], S[1][ind], S[2][ind])
        line = ax.scatter(lst[0], lst[1], lst[2], c = (np.abs(pltlst[i])**2)[ind], s = 5)
        ax.set_xlim(-3, 3)
        ax.set_ylim(-3, 3)
        ax.set_zlim(-3, 3)
        return line,
    return func

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

def overlap(psi1, psi2, space):
    res = np.trapz(np.trapz(np.trapz(np.conj(psi1)*psi2, space[0][:,0][:,0]), space[1][0,:][:,0]), space[2][0,:][0,:])
    return res

print("Error for U3 is {}".format(np.abs(overlap(psi0, psiU3[-1], S)-1)))
print("Error for U72 is {}".format(np.abs(overlap(psi0, psiU72[-1], S)-1)))
print("Error for U11 is {}".format(np.abs(overlap(psi0, psiU11[-1], S)-1)))
print("Error for U7 is {}".format(np.abs(overlap(psi0, psiU7[-1], S)-1)))

pltlst = psiU7
fig = plt.figure()
ax = fig.gca(projection="3d")
ind = np.where(np.abs(pltlst[0])**2 > 0.05)
lst = (S[0][ind], S[1][ind], S[2][ind])
ax.scatter(lst[0], lst[1], lst[2], c = (np.abs(pltlst[0])**2)[ind], s = 5)
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
ax.set_zlim(-3, 3)
ani = animation.FuncAnimation(fig, animate(pltlst, ax), frames=steps, interval = 1, blit = False)

plt.show()
