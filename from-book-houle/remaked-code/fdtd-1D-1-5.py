"""
fdtd-1D-1-5.py
Simulation of a sinusoid wave hitting a lossy dielectric
"""
import numpy as np
from math import pi, sin, exp
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

# Work area width.
width = 200

# Data.
# E-field x-component.
ex = np.zeros(width)

# H-field y-component.
hy = np.zeros(width)

x = np.arange(0, width)

# Pulse parameters
source_position = int(width * 0.1)
t0 = 50
spread = 12

boundary_low = [0, 0]
boundary_high = [0, 0]

# Coefficient corresponds Courant Condition.
# The Courant-Friedrichs-Lewy factor.
cfl_factor = 0.5

c0 = 3e8

# Frequency in MHz
freq_in = 700e6

lambda0 = c0 / freq_in
print(lambda0/10)

# Cell size
ddx = 0.01
ddx = lambda0 / 25
# Time step size
dt = ddx / (2*c0)



# Create Dielectric Profile
epsilon1 = 1
epsilon2 = 4
epsz = 8.854e-12
sigma = 0.04
# sigma = 0


ca = np.empty(width)
ca.fill(epsilon1)

cb = np.empty(width)
cb.fill(epsilon1 * cfl_factor)

cb_start = (int)(width*0.35)
cb_end = (int)(width*0.4)

eaf = dt * sigma / (2 * epsz * epsilon2)
# ca[cb_start:cb_end] = (1 - eaf ) / (1 + eaf)
# cb[cb_start:cb_end] = cfl_factor / (epsilon2 * (1 + eaf))

ca[cb_start:] = (1 - eaf ) / (1 + eaf)
cb[cb_start:] = cfl_factor / (epsilon2 * (1 + eaf))

# ca[cb_start+30:cb_end+30] = (1 - eaf ) / (1 + eaf)
# cb[cb_start+30:cb_end+30] = cfl_factor / (epsilon2 * (1 + eaf))

max_iterations = 10000

# Plot configuration.
fig, (ax1, ax2) = plt.subplots(2)

x_begin = 0
x_end = width
y_begin = -2.0
y_end = 2.0

font_size = 14
plt.rcParams['font.size'] = font_size

ax1.set_xlim(x_begin, x_end)
ax1.set_ylim(y_begin, y_end)
ax1.set_ylabel('E$_x$', fontsize=font_size)

ax2.set_xlim(x_begin, x_end)
ax2.set_ylim(y_begin, y_end)
ax2.set_ylabel('H$_y$', fontsize=font_size)

line1, = ax1.plot([], [], lw=2)
line2, = ax2.plot([], [], lw=2)

line1_profile, = ax1.plot((cfl_factor / cb - 1) / 3, 'k--', linewidth=0.5)
line2_profile, = ax2.plot((cfl_factor / cb - 1) / 3, 'k--', linewidth=0.5)


def init():
    line1.set_data(x, ex)
    line2.set_data(x, ex)
    return line1, line2


def update(iteration_step):

    # Calculate the Ex field.
    for k in range(1, width):
        ex[k] = ca[k] * ex[k] + cb[k] * (hy[k - 1] - hy[k])



    # Electromagnetic "hard" source.
    # Put a Gaussian pulse in the middle.
    pulse = exp(-cfl_factor * ((t0 - iteration_step) / spread) ** 2)* np.sin(freq_in * iteration_step * dt)
    # pulse = exp(-cfl_factor * ((t0 - iteration_step) / spread) ** 2)

    # Put a Sinusoidal "soft"? source.
    # pulse = sin(2*pi*freq_in*dt*iteration_step)

    # Two sources.
    ex[source_position] += pulse
    # ex[int(width * 0.1)] = pulse + ex[int(width * 0.1)]
    # ex[source_position + source_dist] = pulse



    # Absorbing Boundary Conditions
    ex[0] = boundary_low.pop(0)
    boundary_low.append(ex[1])
    ex[width - 1] = boundary_high.pop(0)
    boundary_high.append(ex[width - 2])


    # Calculate the Hy field.
    for k in range(width - 1):
        hy[k] = hy[k] + cfl_factor * (ex[k] - ex[k + 1])

    # Update plot data.
    line1.set_data(x, ex)
    line2.set_data(x, hy)
    title = ax1.text(width / 2, 1.45, 'Iteration step = {}'.format(iteration_step), horizontalalignment='center')

    line1_eps_text = ax1.text(width*0.8, 0.5, 'Eps = {}'.format(epsilon2))
    line2_eps_text = ax2.text(width*0.8, 0.5, 'Eps = {}'.format(epsilon2))

    return line1, line1_profile, line2, line2_profile, title, line1_eps_text, line2_eps_text


anim = FuncAnimation(fig, update, init_func=init,
                     frames=max_iterations, interval=20, blit=True)

plt.show()
