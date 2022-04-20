
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
t = -20:dt:t_max; 

% Speed of light in nm/fs.
c = 299.792458; 

% S = 1 courant factor OLD
% courant factor
cfl_factor=0.98; 

dz = c * dt /  cfl_factor;

% Permeability of free space in (V fs^2/e nm).
mu0=2.013354451e-4; 

% Permittivity of free space in (e / V nm).
eps0=55.26349597e-3; 

% Planck's constant
h = 0.6582119514; 

% Frequency of light in eV (electron-volts).
omega_ev = 1.5d0;

% Frequency of light in (1/fs).
omega = omega_ev / h;

% Width of electric field.
delta = 10.0d0; 

% Electric field is an Gaussian envelop.
source = exp(-(t .^ 2) / (delta^2)) .* cos(omega .* t); 



% Position index of source.
source_position = 50; 

% 1D grid size.
grid_size = 500;

% Space axis
Z = (0:grid_size-1).*dz; 

% Permitivity array.
eps(1:350) = eps0;

% Permitivity array.
eps(350:grid_size) = eps0*1.5; 

% Permeability array.
mu(1:grid_size) = mu0; 

% Conductivity array.
sigma(1:grid_size) = 0;

% Width of PML layer.
pml_width = 55; 

% Polynomial order for grading sigma array (pp 292, Taflove).
m = 3;

neta = sqrt(mu0/eps0);

% required reflectivity
R=1e-8; 

sigma_max =- (m+1) * log(R) / (2 * neta * pml_width * dz);

Pright=((1:pml_width + 1) ./ pml_width) .^ m * sigma_max;

% Lossy conductivity profile.
sigma(grid_size - pml_width:grid_size) = Pright; 
sigma(1:pml_width + 1) = fliplr(Pright);

% Eq 7.8 Taflove, pp 275
sigma_star(1:grid_size) = sigma .* mu0 ./ eps0; 

% PML constants.
A=((mu - 0.5 * dt * sigma_star) ./ (mu + 0.5 * dt * sigma_star)); 
B=(dt / dz) ./ (mu + 0.5 * dt * sigma_star);                          
C=((eps - 0.5 * dt * sigma)./(eps + 0.5 * dt * sigma)); 
D=(dt / dz) ./ (eps + 0.5 * dt * sigma);                     

% Electromagnetic field projections in space array.
Hy(1:grid_size) = 0.0; 
Ex(1:grid_size) = 0.0;

% Plot configuration.
fh = figure(1);
set(fh, 'Color', 'white'); 

% Time loop.
for time_step = 1:length(t)
        % Insert source in certain space grid.
        Ex(source_position) = Ex(source_position) + source(time_step);

        Hy(1:grid_size-1) = A(1:grid_size - 1) .* Hy(1:grid_size - 1) - B(1:grid_size - 1) .* (Ex(2:grid_size) - Ex(1:grid_size - 1));
        Ex(2:grid_size-1) = C(2:grid_size-1) .* Ex(2:grid_size - 1) - D(2:grid_size - 1) .* (Hy(2:grid_size - 1) - Hy(1:grid_size - 2));
        Ex(grid_size) = Ex(grid_size-1);

        % Draw plot.
        figure(1)
        plot(Z,Ex)
        xlabel('Z in nm','fontSize',14);
        ylabel('E_x','fontSize',14);  
        title('1D FDTD with PML','fontSize',14);              
        axis([0 grid_size*dz -1 1]);
        getframe();
end