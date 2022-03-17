"""
    1D fdtd simulation.
    source link: https://www.youtube.com/watch?v=S-6Z8N-30AU
"""

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.animation import FuncAnimation

# Constants
eps0 = 4.85418e-12
mu0 = 1.25664e-6
mu0 = 1.25664e-10

# light speed
c0 = 1 / np.sqrt(eps0 * mu0)

imp0 = np.sqrt(mu0 / eps0)

jmax = 500
jsource = 230
nmax = 2000

Ex = np.zeros(jmax)
Hz = np.zeros(jmax)
x = np.arange(0, jmax)

Ex_prev = np.zeros(jmax)
Hz_prev = np.zeros(jmax)

# Wave lenght in meters.
lambda_min = 350e-9

# Grid steps.
dx = lambda_min / 20
dt = dx / c0

eps = np.ones(jmax) * eps0
eps[200:300] = 4 * eps0
eps[50:100] = 4 * eps0

material_proof = eps > eps0


def source_function(t):
    lambda_0 = 500e-9

    # Frequency
    w0 = 2 * np.pi * c0 / lambda_0
    tau = 30
    t0 = tau * 3

    return np.exp(-(t - t0) ** 2 / tau ** 2) * np.sin(w0 * t * dt)


# Plot configuration.
fig, (ax1, ax2) = plt.subplots(2)

x_begin = 0
x_end = jmax
y_begin = -2.0
y_end = 2.0

font_size = 14
plt.rcParams['font.size'] = font_size

ax1.set_xlim(x_begin, x_end)
ax1.set_ylim(y_begin, y_end)
ax1.set_ylabel('E$_x$', fontsize=font_size)

ax2.set_xlim(x_begin, x_end)
ax2.set_ylim(y_begin, y_end)
ax2.set_ylabel('H$_z$', fontsize=font_size)

line1, = ax1.plot([], [], lw=2)
line2, = ax2.plot([], [], lw=2)

line1_profile = ax1.plot(material_proof)
line2_profile = ax2.plot(material_proof)


def init():
    line1.set_data(x, Ex)
    line2.set_data(x, Hz)
    return line1, line2


def update(i):
    # Update magnetic field boundaries.
    Hz[jmax - 1] = Hz_prev[jmax - 2]

    # Update magnetic field boundaries.
    for j in range(jmax - 1):
        Hz[j] = Hz_prev[j] + dt / (dx * mu0) * (Ex[j + 1] - Ex[j])
        Hz_prev[j] = Hz[j]

    # Magnetic field source.
    Hz[jsource - 1] -= source_function(i) / imp0
    Hz_prev[jsource - 1] = Hz[jsource - 1]

    # Update electric field boundaries.
    Ex[0] = Ex_prev[1]

    # Update electric field.
    for j in range(1, jmax):
        Ex[j] = Ex_prev[j] + dt / (dx * eps[j]) * (Hz[j] - Hz[j - 1])
        Ex_prev[j] = Ex[j]

    # Electric field source.
    Ex[jsource - 1] -= source_function(i + 1)
    Ex_prev[jsource] = Ex[jsource]

    # Update plot data.
    line1.set_data(x, Ex)
    line2.set_data(x, Hz)

    return line1, line2


anim = FuncAnimation(fig, update, init_func=init,
                     frames=nmax, interval=10, blit=True)

plt.show()
