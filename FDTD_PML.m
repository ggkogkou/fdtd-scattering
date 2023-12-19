clear; close; clc;

% Time subdivision to periods
T0 = 25; 

% Maximum Iterations Number
max_steps = 12*T0;

% Speed of light
c = 3e8;

% Magnetic Permeability and Electrical Permitivity of free space
m0 = 4 * pi * 1e-7; % henrys / meter
e0 = 8.854 * 1e-12; % farads / meter
sigma_e = 0;
sigma_m = 0;

% Scaling parameter
scale = 10;

% Grid parameters
cells_x = 5*scale;
cells_y = 10*scale;

% Wave excitation
freq = 10e9; % Hertz
amplitude = 4; % 
source = zeros(1, max_steps); % keep the values of source

% Wavelength
wavelength = c / freq;

% Space and time lattice increment steps
dx = wavelength/10;
dt = dx / (2*c);

% Hard source location
source_x = ceil(0.5*cells_x);
source_y = ceil(0.5*cells_y);

% Simulate source
% Source active time T=10*To, zeros afterwards
for n=1 : 10*T0
    time_n = (n-1) * dt;
    source(n) = amplitude * sin(2 * pi * freq * time_n);
end

% Field arrays
Ez = zeros(cells_x+1, cells_y+1);
Hx = zeros(cells_x+1, cells_y);
Hy = zeros(cells_x, cells_y+1);

% Add metal cylinder
cylinder_diameter = 2*scale;
cylinder_radius = cylinder_diameter / 2;
cylinder_radius_squared = cylinder_radius^2;
cylinder_x = ceil(0.5*cells_x + 3*scale);
cylinder_y = ceil(0.5*cells_y);

cylinder_properties = [3.4, 1, 1.2, 0]; % [er, mr, σe, σm]

% Updating coefficients
Ca = zeros(cells_x+1, cells_y+1);
Cb = zeros(cells_x+1, cells_y+1);
Da_Hx = zeros(cells_x+1, cells_y);
Db_Hx = zeros(cells_x+1, cells_y);
Da_Hy = zeros(cells_x, cells_y+1);
Db_Hy = zeros(cells_x, cells_y+1);

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
            Cb(i, j) = dt/(cylinder_properties(1) * e0 * dx * (1 + tmp));
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

% Define arrays to observe field's variation of two points of the grid
Ez_p1 = zeros(1, max_steps); % p1
Ez_p2 = zeros(1, max_steps); % p2

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
    Ez(source_x, source_y) = source(t);

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

    % Update the values od the two test points p1, p2
    Ez_p1(t) = Ez(scale, scale);
    Ez_p2(t) = Ez(scale, ceil(cells_y/2));

    % Plot
    timestep=int2str(t);
    
    subplot(2,2,1),imagesc(Ez');
    shading flat;
    clim([-1 1]);
    axis([1 cells_x 1 cells_y]);
    colorbar;
    axis image; axis xy
    axis off;
    title(['Ez at time step = ',timestep]);

    % Draw circle around cylinder
    hold on;
    viscircles([cylinder_x, cylinder_y], cylinder_radius, 'EdgeColor', 'm', 'LineWidth', 1);
    hold off;

    Ez_tot(:,:) = Ezx_PML(:,:) + Ezy_PML(:,:);    

    subplot(2,2,2),imagesc(Ez_tot');
    shading flat;
    clim([-2 2]);
    axis([1 cells_x 1 cells_y]);
    colorbar;
    axis image; axis xy
    axis off;
    title(['Ez total inside PML at time step = ',timestep]);

    subplot(2,2,3),imagesc(Hx');
    shading flat;
    clim([-2/377 2/377]);
    axis([1 cells_x+1 1 cells_y]);
    colorbar;
    axis image; axis xy
    axis off;
    title(['Magnetic Field Hx at time step = ',timestep]);

    % Draw circle around cylinder
    hold on;
    viscircles([cylinder_x, cylinder_y], cylinder_radius, 'EdgeColor', 'm', 'LineWidth', 1);
    hold off;
    
    subplot(2,2,4),imagesc(Hy');
    shading flat;
    clim([-2/377 2/377]);
    axis([1 cells_x 1 cells_y+1]);
    colorbar;
    axis image; axis xy
    axis off;
    title(['Magnetic Field Hy at time step = ',timestep]);

    % Draw circle around cylinder
    hold on;
    viscircles([cylinder_x, cylinder_y], cylinder_radius, 'EdgeColor', 'm', 'LineWidth', 1);
    hold off;

    % Save Ez at specific time steps
    %if t == 3*T0 || t == 10*T0 || t == 12*T0
    %    pause;
    %end

    pause(0.01)

end

% Plot the evolution of the field at p1 and p2
figure;
plot(1:max_steps, Ez_p1, 'b-', 1:max_steps, Ez_p2, 'r', 'LineWidth', 2);

title('Plot of Ez_p_1 and Ez_p_2');
xlabel('Steps');
ylabel('Ez values');
legend('Ez_p_1', 'Ez_p_2');
grid on;


