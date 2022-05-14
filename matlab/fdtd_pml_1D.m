% 1D FDTD method.
% With PML Boundary conditions.
% All units in eV,fs, nm.

%--------------------------------------------------------------------%
clc
close all
clear all
%--------------------------------------------------------------------%


% Speed of light in nm/fs.
light_spd = 299.792458;

% S = 1 courant factor OLD
% Courant factor.
cfl_factor = 0.99;

%dz = light_spd * dt /  cfl_factor;
dx = 40;

% Time step in femto seconds.
dt = dx * cfl_factor / light_spd;

% Maximum time limit.
t_max = 300;
t_delay = 10;


% Permeability of free space in (V fs^2/e nm).
mu0 = 2.013354451e-4;

% Permittivity of free space in (e / V nm).
eps0 = 55.26349597e-3;

% Planck's constant.
h = 0.6582119514;

% Frequency of light in eV (electron-volts).
omega_ev = 3.5;

% Frequency of light in (1/fs).
omega = omega_ev / h;

% Width of electric beam.
tau = 4.0;

% Position index of source.
src_position = 75;

% 1D grid size.
grid_size = 500;

% Space axis
X = (0:grid_size - 1) .* dx;

% Permitivity array.
eps(1:349) = eps0;

% Permitivity array.
eps(350:grid_size) = eps0 * 2.5;

% Permeability array.
mu(1:grid_size) = mu0;

% Conductivity array.
sigma(1:grid_size) = 0;

% Width of PML layer.
pml_width = 55;

% Optimal polynomial order for grading sigma array (pp 293, Taflove).
grading_order = 3;

% Impedance.
eta = sqrt(mu0 / eps0);

% Required reflection factor.
R = 1e-8;

% Taflove, pp 293, Eq 7.62.
sigma_max = -(grading_order+1) * log(R) / (2 * eta * pml_width * dx);

% Taflove, pp 292, Eq 7.60a.
for i = 1:pml_width
    sigma_in_pml = ((i / pml_width) .^ grading_order) * sigma_max;

    % Lossy electric conductivity profile.
    sigma(grid_size-pml_width+i) = sigma_in_pml;
    sigma(pml_width-i+1) = sigma_in_pml;
end


% Eq 7.8 Taflove, pp 275
% Magnetic conductivity loss.
sigma_star(1:grid_size) = sigma .* mu0 ./ eps0;

% PML constants.
A = ((mu - 0.5 * dt * sigma_star) ./ (mu + 0.5 * dt * sigma_star));
B = (dt / dx) ./ (mu + 0.5 * dt * sigma_star);
C = ((eps - 0.5 * dt * sigma)./(eps + 0.5 * dt * sigma));
D = (dt / dx) ./ (eps + 0.5 * dt * sigma);

% Electromagnetic field projections in space array.
hy(1:grid_size) = 0.0;
ex(1:grid_size) = 0.0;

% Plot configuration.
fh = figure(1);
set(fh, 'Color', 'white');

% Time loop.
for time_step = 1:t_max
        % Insert source in certain space grid.
        % Electric field is an Gaussian envelop.
        t = time_step * dt;
        source = exp(-(((t_delay - t) / tau) ^ 2)) * cos(omega * t);
        ex(src_position) = ex(src_position) + source;


        % Hy field calculation.
        for i = 1:grid_size - 1
            hy(i) = A(i) * hy(i) - B(i) * (ex(i+1) - ex(i));
        end


        % Ex field calculation.
        % i= 1 = pec, so don't evaluate.
        for i = 2:grid_size - 1
            ex(i) = C(i) * ex(i) - D(i) * (hy(i) - hy(i-1));
        end

        %ex(grid_size) = ex(grid_size - 1);


        % Draw plot.
        figure(1)
        plot(X, ex)
        xlabel('Z in nm','fontSize', 14);
        ylabel('E_x','fontSize', 14);

        time_step_str = int2str(time_step);
        title(['1D FDTD with PML. Time step = ', time_step_str],'fontSize', 14);
        axis([0 grid_size*dx -1 1]);
        getframe();
end