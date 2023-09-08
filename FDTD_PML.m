clear; close; clc;

% Maximum Iterations Number
max_steps = 300;

% Speed of light
c = 3e8;

% Magnetic Permeability and Electrical Permitivity of free space
m0 = 4 * pi * 1e-7; % henrys / meter
e0 = 8.854 * 1e-12; % farads / meter
sigma_e = 0;
sigma_m = 0;

% Grid parameters
cells_x = 100;
cells_y = 100;

% Space and time lattice increment steps
dx = 3e-3;
dt = dx / (2*c);

% Hard source location
hard_source_x = ceil(cells_x/10);
hard_source_y = ceil(cells_y/2);

% Wave excitation
freq = 5e9; % Hertz
amplitude = 2; % 
hard_source = zeros(1, max_steps); % keep the values of source

for t=1 : max_steps
    time_t = (t-1) * dt;
    hard_source(t) = amplitude * sin(2 * pi * freq * time_t);
end

% Field arrays
Ez = zeros(cells_x+1, cells_y+1);
Hx = zeros(cells_x, cells_y+1);
Hy = zeros(cells_x+1, cells_y);

% Add metal cylinder
cylinder_diameter = 20;
cylinder_radius = cylinder_diameter / 2;
cylinder_radius_squared = cylinder_radius^2;
cylinder_x = ceil(0.8 * cells_x);
cylinder_y = ceil(0.5 * cells_y);

cylinder_properties = [1, 1, 1e7, 0]; % [er, mr, σe, σm]

% Updating coefficients
Ca = zeros(cells_x+1, cells_y+1);
Cb = zeros(cells_x+1, cells_y+1);
Da_Hx = zeros(cells_x, cells_y+1);
Db_Hx = zeros(cells_x, cells_y+1);
Da_Hy = zeros(cells_x+1, cells_y);
Db_Hy = zeros(cells_x+1, cells_y);

tmp = sigma_e*dt/(2*e0);
Ca(:, :) = (1-tmp)/(1+tmp);
Cb(:, :) = dt/e0/dx/(1+tmp);

tmp = sigma_m*dt/(2*m0);
Da_Hx(:, :) = (1-tmp)/(1+tmp);
Db_Hx(:, :) = dt/m0/dx/(1+tmp);
Da_Hy(:, :) = (1-tmp)/(1+tmp);
Db_Hy(:, :) = dt/m0/dx/(1+tmp);

% Add the cylinder - optimized grid iteration
% Search only inside the cylinder area, not the whole grid
x1 = max(1, floor(cylinder_x - cylinder_radius));
x2 = min(cells_x, ceil(cylinder_x + cylinder_radius));
y1 = max(1, floor(cylinder_y - cylinder_radius));
y2 = min(cells_y, ceil(cylinder_y + cylinder_radius));

% Loop through the relevant grid cells
for i=x1 : x2
    for j=y1 : y2
        % Calculate squared distance to the cylinder center
        dist_squared = (i - cylinder_x)^2 + (j - cylinder_y)^2;
        
        % Check if the cell is inside the cylinder
        if dist_squared <= cylinder_radius_squared
            % Assign values for electric field
            tmp = cylinder_properties(3)*dt/(2*cylinder_properties(1)*e0);
            Ca(i, j) = (1-tmp)/(1+tmp);
            Cb(i, j) = dt/cylinder_properties(1)/dx/(1+tmp);
        end
    end
end

% PML parameters + arrays + conductivities (polynomial grading)
pml_cells = 8;
pml_R0 = 1e-7;
pml_order = 2;

Ezx_PML = zeros(pml_cells, cells_y+1);
Ezy_PML = zeros(pml_cells, cells_y+1);
Hx_PML = zeros(pml_cells, cells_y);
Hy_PML = zeros(pml_cells, cells_y+1);

Ca_PML = zeros(pml_cells, cells_y+1);
Cb_PML = zeros(pml_cells, cells_y+1);

Da_Hx_PML = zeros(pml_cells, cells_y);
Db_Hx_PML = zeros(pml_cells, cells_y);
Da_Hy_PML = zeros(pml_cells, cells_y+1);
Db_Hy_PML = zeros(pml_cells, cells_y+1);

sigma_max = -(pml_order+1)*e0*c*log(pml_R0)/(2*pml_cells*dx);

% Find coefficients Ca, Cb, Da, Db for the PML region
for x=2 : pml_cells
    x1 = pml_cells - x + 0.5;
    x2 = pml_cells - x + 1.5;
    sigma_x_e = sigma_max / ((pml_order+1)*pml_cells^pml_order) * (x2^(pml_order+1)-x1^(pml_order+1));
    sigma_x_m = sigma_x_e * m0 / e0;

    Ca_PML(x, 2:cells_y) = exp(-sigma_x_e*dt/e0);
    Cb_PML(x, 2:cells_y) = (1-Ca_PML(x, 2:cells_y))/(sigma_x_e*dx);

    Da_Hx_PML(x, 1:cells_y) = exp(- sigma_x_m*dt/m0);
    Db_Hx_PML(x, 1:cells_y) = (1-Da_Hx_PML(x, 1:cells_y))/(sigma_x_m*dx);

    x1 = pml_cells - x;
    x2 = pml_cells - x + 1;
    sigma_x_e = sigma_max / ((pml_order+1)*pml_cells^pml_order) * (x2^(pml_order+1)-x1^(pml_order+1));
    sigma_x_m = sigma_x_e * m0 / e0;

    Da_Hy_PML(x, 2:cells_y) = exp(-sigma_x_m*dt/m0);
    Db_Hy_PML(x, 2:cells_y) = (1-Da_Hy_PML(x, 2:cells_y))/(sigma_x_m*dx);
end

x1 = pml_cells-1;
x2 = pml_cells;
sigma_x_e = sigma_max / ((pml_order+1)*pml_cells^pml_order) * (x2^(pml_order+1)-x1^(pml_order+1));
sigma_x_m = sigma_x_e * m0 / e0;

Da_Hy_PML(1, 2:cells_y) = exp(-sigma_x_m*dt/m0);
Db_Hy_PML(1, 2:cells_y) = (1-Da_Hy_PML(1, 2:cells_y))/(sigma_x_m*dx);

% Main FDTD Loop
for t=1 : max_steps

    % Update Ez (Vacuum)
    for x=2 : cells_x
        for y=2 : cells_y
            Ez(x, y) = Ca(x, y)*Ez(x, y) + Cb(x, y)*(Hy(x, y)-Hy(x-1,y)+Hx(x, y-1)-Hx(x, y));
        end
    end

    % Update boundary
    for y=2 : cells_y
        Ez(1, y) = Ca(1, y)*Ez(1, y) + Cb(1, y)*(Hy(1, y)-Hy_PML(pml_cells, y)+Hx(1, y-1)-Hx(1, y));
    end

    % Source excitation
    Ez(hard_source_x, hard_source_y) = hard_source(t);

    % Update Ezx, Ezy (PML)
    for x=2 : pml_cells
        for y=2 : cells_y
            Ezx_PML(x, y) = Ca_PML(x, y)*Ezx_PML(x, y) + Cb_PML(x, y)*(Hy_PML(x, y)-Hy_PML(x-1, y));
            Ezy_PML(x, y) = Ca_PML(x, y)*Ezy_PML(x, y) + Cb_PML(x, y)*(Hx_PML(x, y-1)-Hx_PML(x, y));
        end
    end

    % Update Hx (Vacuum)
    for x=1 : cells_x
        for y=1 : cells_y
            Hx(x, y) = Da_Hx(x, y)*Hx(x, y) + Db_Hx(x, y)*(Ez(x, y)-Ez(x,y+1));
        end
    end

    % Update Hx (PML)
    for x=2 : pml_cells
        for y=1 : cells_y
            Hx_PML(x, y) = Da_Hx_PML(x, y)*Hx_PML(x, y) + Db_Hx_PML(x, y)*(Ezx_PML(x, y)+Ezy_PML(x, y)-Ezx_PML(x, y+1)-Ezy_PML(x, y+1));
        end
    end

    % Update Hy (Vacuum)
    for x=1 : cells_x
        for y=2 : cells_y
            Hy(x, y) = Da_Hy(x, y)*Hy(x, y) + Db_Hy(x, y)*(Ez(x+1, y)-Ez(x, y));
        end
    end

    % Update Hy (PML)
    for x=1 : pml_cells-1
        for y=2 : cells_y
            Hy_PML(x, y) = Da_Hy_PML(x, y)*Hy_PML(x, y) + Db_Hy_PML(x, y)*(Ezx_PML(x+1, y)+Ezy_PML(x+1, y)-Ezx_PML(x, y)-Ezy_PML(x, y));
        end
    end

    for y=2 : cells_y
        Hy_PML(pml_cells, y) = Da_Hy_PML(pml_cells, y)*Hy_PML(pml_cells, y) + Db_Hy_PML(pml_cells, y)*(Ez(1, y)-Ezx_PML(pml_cells, y)-Ezy_PML(pml_cells, y));
    end

    % Plot
    timestep=int2str(t);
    
    subplot(1,2,1),imagesc(Ez');
    shading flat;
    clim([-1 1]);
    axis([1 cells_x 1 cells_y]);
    colorbar;
    axis image; axis xy
    axis off;
    title(['Ez at time step = ',timestep]);

    subplot(1,2,2),imagesc(Ezx_PML');
    shading flat;
    clim([-1 1]);
    axis([1 cells_x 1 cells_y]);
    colorbar;
    axis image; axis xy
    axis off;
    title(['Ez at time step = ',timestep]);

    pause(0.01)

end

