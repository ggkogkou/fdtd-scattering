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
cells_x = 10*scale;
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
Hx = zeros(cells_x, cells_y+1);
Hy = zeros(cells_x+1, cells_y);

% Add metal cylinder
cylinder_diameter = 2*scale;
cylinder_radius = cylinder_diameter / 2;
cylinder_radius_squared = cylinder_radius^2;
cylinder_x = ceil(0.5*cells_x + 3*scale);
cylinder_y = ceil(0.5*cells_y);

cylinder_properties = [3.4, 1, 1.2, 0]; % [er, mr, σe, σm]

% Updating coefficients
Ca = zeros(cells_x, cells_y);
Cb = zeros(cells_x, cells_y);
Da = zeros(cells_x, cells_y);
Db = zeros(cells_x, cells_y);

tmp = sigma_e*dt/(2*e0);
Ca(:, :) = (1-tmp)/(1+tmp);
Cb(:, :) = dt/e0/dx/(1+tmp);

tmp = sigma_m*dt/(2*m0);
Da(:, :) = (1-tmp)/(1+tmp);
Db(:, :) = dt/m0/dx/(1+tmp);

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
            Cb(i, j) = dt/cylinder_properties(1)/e0/dx/(1+tmp);
        end
    end
end

% Main FD-TD Loop
for t=1 : max_steps

    % Store Ez(2, :) values for Mur's ABC
    tmp_mur = Ez(2, :);

    % Update Ez
    for x=2 : cells_x
        for y=2 : cells_y
            Ez(x, y) = Ca(x, y)*Ez(x, y) + Cb(x, y)*(Hy(x, y)-Hy(x-1,y)+Hx(x, y-1)-Hx(x, y));
        end
    end

    % Update left side boundaries using Mur's 1st Order ABC
    Ez(1, :) = tmp_mur(1, :) - (dx-c*dt)/(dx+c*dt)*(Ez(2, :)-Ez(1, :));

    % Source excitation
    Ez(source_x, source_y) = source(t);

    % Update Hx
    for x=2 : cells_x
        for y=1 : cells_y
            Hx(x, y) = Da(x, y)*Hx(x, y) + Db(x, y)*(Ez(x, y)-Ez(x,y+1));
        end
    end

    % Update Hy
    for x=1 : cells_x
        for y=2 :cells_y
            Hy(x, y) = Da(x, y)*Hy(x, y) + Db(x, y)*(Ez(x+1, y)-Ez(x, y));
        end
    end

    % Plot
    timestep=int2str(t);
    
    subplot(1,1,1),imagesc(Ez');
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

    % Save Ez at specific time steps
    %if t == 3*T0 || t == 5*T0 || t == 8*T0
    %   pause;
    %end
    
    pause(0.01)
    
end

