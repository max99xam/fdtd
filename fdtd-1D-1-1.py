"""
fdtd-1D-1-1.py
Simulation in free space.
"""
import numpy as np
from math import exp
from matplotlib import pyplot as plt

width = 250
ex = np.zeros(width)
hy = np.zeros(width)

# Pulse parameters
kc = int(width / 2)
t0 = 40
spread = 12

nsteps = 100

# Main FDTD Loop
for time_step in range(1, nsteps + 1):
    # Calculate the Ex field
    for k in range(1, width):
        ex[k] = ex[k] + 0.5 * (hy[k - 1] - hy[k])

    # Electromagnetic "hard" source.
    # Put a Gaussian pulse in the middle
    pulse = exp(-0.5 * ((t0 - time_step) / spread) ** 2)

    source_dist = 50

    # Two sources.
    ex[kc-source_dist] = pulse
    ex[kc+source_dist] = pulse

    # Calculate the Hy field
    for k in range(width - 1):
        hy[k] = hy[k] + 0.5 * (ex[k] - ex[k + 1])

# Plot Ex
plt.rcParams['font.size'] = 12
plt.figure(figsize=(8, 3.5))
plt.subplot(211)
plt.plot(ex, color='red', linewidth=1)
plt.ylabel('E$_x$', fontsize='14')
plt.xticks(np.arange(0, width+1, step=20))
plt.xlim(0, width)
plt.yticks(np.arange(-1, 1.2, step=1))
plt.ylim(-1.2, 1.2)
plt.text(width/2, 1.5, 'Time step = {}'.format(time_step), horizontalalignment='center')

# Plot Hy
plt.subplot(212)
plt.plot(hy, color='k', linewidth=1)
plt.ylabel('H$_y$', fontsize='14')
plt.xlabel('FDTD cells')
plt.xticks(np.arange(0, width+1, step=20))
plt.xlim(0, width)
plt.yticks(np.arange(-1, 1.2, step=1))
plt.ylim(-1.2, 1.2)

plt.subplots_adjust(bottom=0.2, hspace=0.45)
plt.show()
