% 2D FDTD with PML.

% Based on code from book:
% Electromagnetic simulation using the fdtd method with python.
% Chapter 3.

% Grid sizes.
pml_width = 40;
rows = 200 + pml_width * 2;
cols = 200 + pml_width * 2;

% Set source position
src_row = rows / 2;
src_col = cols / 2;

% Light speed.
light_spd = 2.99792458e8;

% Permeability of free space.
mu0 = 4.0 * pi * 1.0e-7;

% Permittivity of free space.
epsz = 8.8e-12;

% Permittivity.
epsilon0 = 1;
epsilon1 = 1;
epsilon2 = 3.2;


% Conductivity.
sigma0 = 5000;
sigma1 = 0.001;
sigma2 = 0.001;

dz = zeros(rows,cols);
ez = zeros(rows,cols);
hx = zeros(rows,cols);
hy = zeros(rows,cols);

gaz = ones(rows,cols);

% Max simulation time.
max_time = 400;

% Space grid step.
ddx = 1.0e-3;
ddy = ddx;

% Courant factor.
cfl_factor = 0.98;

% Time step corressponding Courant factor.
dt = cfl_factor / (light_spd * sqrt((1 / ddx)^2 + (1 / ddy)^2));

% Source params.
% Gaussian beam.
t0 = 20;

% Beam width.
tau = 20;

% PML params.
gi2 = ones(rows);
gi3 = ones(rows);
fi1 = zeros(rows);
fi2 = ones(rows);
fi3 = ones(rows);

gj2 = ones(rows);
gj3 = ones(rows);
fj1 = zeros(rows);
fj2 = ones(rows);
fj3 = ones(rows);

% Calculate PML.
for i = 1:pml_width
    xnum = pml_width - i;
    xxn = xnum / pml_width;
    xn = 0.33 * ((xxn)^3);

    gi2(i) = 1.0 / (1.0 + xn);
    gi2(rows - 1 - i) = 1.0 / (1.0 + xn);
    gi3(i) = (1.0 - xn) / (1.0 + xn);
    gi3(rows - i - 1) = (1.0 - xn) / (1.0 + xn);

    gj2(i) = 1.0 / (1.0 + xn);
    gj2(rows - 1 - i) = 1.0 / (1.0 + xn);
    gj3(i) = (1.0 - xn) / (1.0 + xn);
    gj3(rows - i - 1) = (1.0 - xn) / (1.0 + xn);

    xxn = (xnum - 0.5) / pml_width;
    xn = 0.33 * ((xxn)^3);

    fi1(i) = xn;
    fi1(rows - 2 - i) = xn;
    fi2(i) = 1.0 / (1.0 + xn);
    fi2(rows - 2 - i) = 1.0 / (1.0 + xn);
    fi3(i) = (1.0 - xn) / (1.0 - xn);
    fi3(rows - 2 - i) = (1.0 - xn) / (1.0 + xn);

    fj1(i) = xn;
    fj1(rows - 2 - i) = xn;
    fj2(i) = 1.0 / (1.0 + xn);
    fj2(rows - 2 - i) = 1.0 / (1.0 + xn);
    fj3(i) = (1.0 - xn) / (1.0 - xn);
    fj3(rows - 2 - i) = (1.0 - xn) / (1.0 + xn);
end

for i = 1:rows
    for j = 1:cols

        % Medium 1.
        gaz(i,j) = 1 ./ (epsilon1 + (sigma1 * dt) / epsz);

        % Medium 2.
        if i>=90 && i<=130 && j>=90 && j<=130
            gaz(i,j) = 1 ./ (epsilon2 + (sigma2 * dt) / epsz);
        end

        % Dielectric border. Medium 0.
        %if i>=(src_row - 20) && i<=(src_row + 20) ...
        %   && (j == (src_col-10) || j == (src_col+10)) ...
        %   || j>=(src_col - 10) && j<=(src_col + 10) ...
        %   && (i == (src_row-10))
        %        gaz(i,j) = 1 ./ (epsilon0 + (sigma0 * dt) / epsz);
        %end
    end
end


% Main fdtd loop.
% Time layers calculation.
for time_step = 1:max_time

    % Dz field calculation.
    for j = 2:cols
        for i = 2:rows
            dz(i,j) = gi3(i) * gj3(j) * dz(i,j) ...
                    + gi2(i) * gj2(j) * .5 ...
                    * (hy(i,j) - hy(i-1,j) - hx(i,j) + hx(i,j-1));
        end
    end

    % Ez field calculation.
    for j = 2:cols
        for i = 2:rows
            ez(i,j) = gaz(i,j) * dz(i,j);
        end
    end

    % Put Gaussian beam source.
    source = -2.0 * ((time_step - t0) ./ tau) .* exp(-1.0 * ((time_step - t0) ./ tau) .^ 2.0);
    ez(src_row, src_col) = source;

    % Hx field calculation.
    for j = 1:cols-1
        for i = 1:rows-1
            hx(i,j) = fj3(j) * hx(i,j) + fj2(j) * 0.5 * (ez(i,j) - ez(i,j+1));
        end
    end

    % Hy field calculation.
    for j = 1:cols-1
        for i = 1:rows-1
            hy(i,j) = fi3(i) * hy(i,j) + fi2(i) * 0.5 * (ez(i+1,j) - ez(i,j));
        end
    end

    max_ez = max(ez);
    min_ez = min(ez);

    disp("~~~~~~~~~~")
    disp(max_ez)
    disp(min_ez)
    disp("----------")

    % Draw plot.
    time_step_str = int2str(time_step);
    
    min_val = -0.05;
    max_val = 0.05;
    ez_lims = [min_val max_val];
    imagesc(ez, ez_lims);

    colormap(jet);
    colorbar
    title(['Ez at time step = ', time_step_str]);

    pause(0.001);
end

