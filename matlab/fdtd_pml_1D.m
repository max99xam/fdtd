
% 1D FDTD method.
% With PML Boundary conditions.
% All units in eV,fs, nm.

%--------------------------------------------------------------------%
clc
close all
clear all
%--------------------------------------------------------------------%

% Time step in femto seconds.
dt = 0.13; 

% Maximum time limit.
t_max = 1000; 

% Time array.
t = 0:dt:t_max; 
t0 = 20;

% Speed of light in nm/fs.
light_spd = 299.792458; 

% S = 1 courant factor OLD
% Courant factor.
cfl_factor = 0.98;

dz = light_spd * dt /  cfl_factor;

% Permeability of free space in (V fs^2/e nm).
mu0 = 2.013354451e-4;

% Permittivity of free space in (e / V nm).
eps0 = 55.26349597e-3;

% Planck's constant.
h = 0.6582119514;

% Frequency of light in eV (electron-volts).
omega_ev = 1.5d0;

% Frequency of light in (1/fs).
omega = omega_ev / h;

% Width of electric beam.
tau = 10.0d0;

% Electric field is an Gaussian envelop.
source = exp(-(((t0 - t) / tau) .^ 2)) .* cos(omega .* t);

% Position index of source.
src_position = 150;

% 1D grid size.
grid_size = 500;

% Space axis
Z = (0:grid_size - 1) .* dz;

% Permitivity array.
eps(1:350) = eps0;

% Permitivity array.
eps(350:grid_size) = eps0 * 1.5;

% Permeability array.
mu(1:grid_size) = mu0;

% Conductivity array.
sigma(1:grid_size) = 0;

% Width of PML layer.
pml_width = 55;

% Polynomial order for grading sigma array (pp 292, Taflove).
m = 3;

% Impedance.
eta = sqrt(mu0 / eps0);

% Required reflection factor.
R = 1e-8;

% Taflove, pp 292, Eq 7.61.
sigma_max =- (m+1) * log(R) / (2 * eta * pml_width * dz);

% Taflove, pp 292, Eq 7.60a.
Pright = ((1:pml_width + 1) ./ pml_width) .^ m * sigma_max;

% Lossy electric conductivity profile.
sigma(grid_size - pml_width:grid_size) = Pright;
sigma(1:pml_width + 1) = fliplr(Pright);

% Eq 7.8 Taflove, pp 275
% Magnetic conductivity loss.
sigma_star(1:grid_size) = sigma .* mu0 ./ eps0;

% PML constants.
A = ((mu - 0.5 * dt * sigma_star) ./ (mu + 0.5 * dt * sigma_star));
B = (dt / dz) ./ (mu + 0.5 * dt * sigma_star);
C = ((eps - 0.5 * dt * sigma)./(eps + 0.5 * dt * sigma));
D = (dt / dz) ./ (eps + 0.5 * dt * sigma);

% Electromagnetic field projections in space array.
hy(1:grid_size) = 0.0;
ex(1:grid_size) = 0.0;

% Plot configuration.
fh = figure(1);
set(fh, 'Color', 'white');

% Time loop.
for time_step = 1:length(t)
        % Insert source in certain space grid.
        ex(src_position) = ex(src_position) + source(time_step);

        hy(1:grid_size - 1) = A(1:grid_size - 1) .* hy(1:grid_size - 1) - B(1:grid_size - 1) .* (ex(2:grid_size) - ex(1:grid_size - 1));
        ex(2:grid_size - 1) = C(2:grid_size - 1) .* ex(2:grid_size - 1) - D(2:grid_size - 1) .* (hy(2:grid_size - 1) - hy(1:grid_size - 2));
        ex(grid_size) = ex(grid_size - 1);

        % Draw plot.
        figure(1)
        plot(Z, ex)
        xlabel('Z in nm','fontSize', 14);
        ylabel('E_x','fontSize', 14);
        title('1D FDTD with PML','fontSize', 14);
        axis([0 grid_size*dz -1 1]);
        getframe();
end