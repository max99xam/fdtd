'''
1D FDTD method.
With PML Boundary conditions.
All units in eV,fs, nm.
'''
import numpy as np
from math import pi, sin, exp
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


# Speed of light in nm/fs.
light_spd = 299.792458

# S = 1 courant factor OLD
# Courant factor.
cfl_factor = 0.99

#dz = light_spd * dt /  cfl_factor
dx = 40

# Time step in femto seconds.
dt = dx * cfl_factor / light_spd

# Maximum time limit.
t_max = 1300
t_delay = 8


# Permeability of free space in (V fs^2/e nm).
mu0 = 2.013354451e-4

# Permittivity of free space in (e / V nm).
eps0 = 55.26349597e-3

# Planck's constant.
h = 0.6582119514

# Frequency of light in eV (electron-volts).
omega_ev = 1.25

# Frequency of light in (1/fs).
omega = omega_ev / h

# Width of electric beam.
tau = 6.0

# Position index of source.
src_position = 100

# 1D grid size.
grid_size = 500

# Space axis
#X = (0:grid_size - 1) .* dx
x = np.arange(0, grid_size)


# Permitivity array.
eps = np.empty(grid_size)
eps.fill(eps0)
eps[300:grid_size] = eps0*2


# Permeability array.
mu = np.empty(grid_size)
mu.fill(mu0)

# Conductivity array.
sigma = np.empty(grid_size)
sigma.fill(0)

# Width of PML layer.
pml_width = 60

# Optimal polynomial order for grading sigma array (pp 293, Taflove).
grading_order = 3

# Impedance.
eta = np.sqrt(mu0 / eps0)

# Required reflection factor.
R = 1e-8

# Taflove, pp 293, Eq 7.62.
sigma_max = -(grading_order+1) * np.log(R) / (2 * eta * pml_width * dx)

# Taflove, pp 292, Eq 7.60a.
for i in range(pml_width):
    sigma_in_pml = np.power(i / pml_width, grading_order) * sigma_max

    # Lossy electric conductivity profile.
    sigma[grid_size-pml_width+i] = sigma_in_pml
    sigma[pml_width-i-1] = sigma_in_pml


# Eq 7.8 Taflove, pp 275
# Magnetic conductivity loss.
sigma_star = np.empty(grid_size)
A = np.empty(grid_size)
B = np.empty(grid_size)
C = np.empty(grid_size)
D = np.empty(grid_size)

for i in range(grid_size):

    # Eq 7.8 Taflove, pp 275
    # Magnetic conductivity loss.
    sigma_star[i] = sigma[i] * mu0 / eps0

    # PML constants.
    A[i] = ((mu[i] - 0.5 * dt * sigma_star[i]) / (mu[i] + 0.5 * dt * sigma_star[i]))

    B[i] = (dt / dx) / (mu[i] + 0.5 * dt * sigma_star[i])

    C[i] = ((eps[i] - 0.5 * dt * sigma[i]) / (eps[i] + 0.5 * dt * sigma[i]))

    D[i] = (dt / dx) / (eps[i] + 0.5 * dt * sigma[i])


# Electromagnetic field projections in space array.
hy = np.empty(grid_size)
hy.fill(0)

ex = np.empty(grid_size)
ex.fill(0)

# Plot configuration.
fig, (ax1, ax2) = plt.subplots(2)

x_begin = 0
x_end = grid_size
y_begin = -1.0
y_end = 1.0

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


# Time loop.
def update(time_step):
        # Insert source in certain space grid.
        # Electric field is an Gaussian envelop.
        t = time_step * dt
        source = np.exp(-(np.power(t_delay - t / tau, 2))) * np.cos(omega * t)
        ex[src_position] = ex[src_position] + source


        # Hy field calculation.
        for i in range(grid_size - 1):
            hy[i] = A[i] * hy[i] - B[i] * (ex[i+1] - ex[i])


        # Ex field calculation.
        # i= 1 = pec, so don't evaluate.
        for i in range(1,grid_size - 1):
            ex[i] = C[i] * ex[i] - D[i] * (hy[i] - hy[i-1])

        # Update plot data.
        line1.set_data(x, ex)
        line2.set_data(x, hy)
        title = ax1.text(grid_size / 2, 1.45, 'Time step = {}'.format(time_step), horizontalalignment='center')

        # line1_eps_text = ax1.text(width*0.8, 0.5, 'Eps = {}'.format(epsilon2))
        # line2_eps_text = ax2.text(width*0.8, 0.5, 'Eps = {}'.format(epsilon2))

        return line1, line2, title

anim = FuncAnimation(fig, update, init_func=init,
                     frames=t_max, interval=20, blit=True)

plt.show()