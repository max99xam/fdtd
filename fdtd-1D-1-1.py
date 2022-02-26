"""
fdtd-1D-1-1.py
Simulation in free space (without boundary condition).
"""
import numpy as np
from math import exp
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

# Work area width.
width = 150

# Data.
ex = np.zeros(width)
hy = np.zeros(width)
x = np.arange(0, width)

# Pulse parameters
source_position = int(width / 2)
t0 = 40
spread = 12

# Coefficient corresponds Courant Condition.
coeff = 0.5

max_iterations = 1200

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


def init():
    line1.set_data(x, ex)
    line2.set_data(x, ex)
    return line1, line2


def update(iteration_step):
    # Calculate the Ex field.
    for k in range(1, width):
        ex[k] = ex[k] + coeff * (hy[k - 1] - hy[k])

    # Electromagnetic "hard" source.
    # Put a Gaussian pulse in the middle.
    pulse = exp(-coeff * ((t0 - iteration_step) / spread) ** 2)

    source_dist = 20

    # Two sources.
    ex[source_position - source_dist] = pulse
    ex[source_position + source_dist] = pulse

    # Calculate the Hy field.
    for k in range(width - 1):
        hy[k] = hy[k] + coeff * (ex[k] - ex[k + 1])

    # Update plot data.
    line1.set_data(x, ex)
    line2.set_data(x, hy)
    title = ax1.text(width / 2, 1.45, 'Iteration step = {}'.format(iteration_step), horizontalalignment='center')

    return line1, line2, title


anim = FuncAnimation(fig, update, init_func=init,
                     frames=max_iterations, interval=20, blit=True)

plt.show()
