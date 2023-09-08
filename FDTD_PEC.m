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

% Cell array that stores all the different media properties
media = {[1, 1, sigma_e, sigma_m]};

% Grid parameters
total_cells_x = 90;
total_cells_y = 40;

% Space and time lattice increment steps
dx = 3e-3;
dt = dx / (2*c);

% Hard source location
hard_source_x = ceil(total_cells_x/3);
hard_source_y = ceil(total_cells_y/2);

% Wave excitation
freq = 5e9; % Hertz
amplitude = 2; % 
hard_source = zeros(1, max_steps); % keep the values of source

for n=1 : max_steps
    time_n = (n-1) * dt;
    hard_source(n) = amplitude * sin(2 * pi * freq * time_n);
end

% Field arrays
Ez = zeros(total_cells_x+1, total_cells_y+1);
Hx = zeros(total_cells_x, total_cells_y+1);
Hy = zeros(total_cells_x+1, total_cells_y);

% Add metal cylinder
cylinder_diameter = 20;
cylinder_radius = cylinder_diameter / 2;
cylinder_radius_squared = cylinder_radius^2;
cylinder_x = ceil(0.8 * total_cells_x);
cylinder_y = ceil(0.5 * total_cells_y);

cylinder_properties = [1, 1, 1e7, 0]; % [er, mr, σe, σm]

% Add cylinder as second media
media{2} = cylinder_properties;

% Updating coefficients
number_of_media = size(media, 2);
Ca = zeros(total_cells_x, total_cells_y);
Cb = zeros(total_cells_x, total_cells_y);
Da = zeros(total_cells_x, total_cells_y);
Db = zeros(total_cells_x, total_cells_y);

tmp = sigma_e*dt/(2*e0);
Ca(:, :) = (1-tmp)/(1+tmp);
Cb(:, :) = dt/e0/dx/(1+tmp);

tmp = sigma_m*dt/(2*m0);
Da(:, :) = (1-tmp)/(1+tmp);
Db(:, :) = dt/m0/dx/(1+tmp);

% Add the cylinder - optimized grid iteration
% Search only inside the cylinder area, not the whole grid
x1 = max(1, floor(cylinder_x - cylinder_radius));
x2 = min(total_cells_x, ceil(cylinder_x + cylinder_radius));
y1 = max(1, floor(cylinder_y - cylinder_radius));
y2 = min(total_cells_y, ceil(cylinder_y + cylinder_radius));

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

% Main FD-TD Loop
for t=1 : max_steps

    % Update Ez
    for x=2 : total_cells_x
        for y=2 : total_cells_y
            Ez(x, y) = Ca(x, y)*Ez(x, y) + Cb(x, y)*(Hy(x, y)-Hy(x-1,y)+Hx(x, y-1)-Hx(x, y));
        end
    end

    % Source excitation
    Ez(hard_source_x, hard_source_y) = hard_source(t);

    % Update Hx
    for x=2 : total_cells_x
        for y=1 : total_cells_y
            Hx(x, y) = Da(x, y)*Hx(x, y) + Db(x, y)*(Ez(x, y)-Ez(x,y+1));
        end
    end

    % Update Hy
    for x=1 : total_cells_x
        for y=2 :total_cells_y
            Hy(x, y) = Da(x, y)*Hy(x, y) + Db(x, y)*(Ez(x+1, y)-Ez(x, y));
        end
    end

    % Plot
    timestep=int2str(t);
    
    subplot(1,1,1),imagesc(Ez');
    shading flat;
    clim([-1 1]);
    axis([1 total_cells_x 1 total_cells_y]);
    colorbar;
    axis image; axis xy
    axis off;
    title(['Ez at time step = ',timestep]);
    pause(0.01)

end

