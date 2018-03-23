function plotting_output = pic_sub_driver(problem, Grid_Size, Remapping)


%Grid_Size should start around 8 and end being doubled at 32 
%Grid_Size = 8;

%Parameters from papers for distribution and size of grid
k = 0.5;
alpha = 0.01;
vmax = 10;
L = 2*pi/k;

dt = 1/(1*Grid_Size);
T = 30;
NumCells = 2*Grid_Size;
Numx = 4*Grid_Size;
Numv = 8*Grid_Size;



%Spatial Grid
spatial_grid_init = transpose(linspace(0,L,NumCells+1));
grid.poisson_grid = spatial_grid_init(1:end-1, :);

%Grids along x and v axis with particles at even placements
particle_spatial_grid_init = transpose(linspace(0,L,Numx+1));
grid.x_grid = particle_spatial_grid_init(1:end-1, :);
particle_velocity_grid_init = transpose(linspace(-vmax, vmax,Numv+1));
grid.v_grid = particle_velocity_grid_init(1:end-1, :);


hx = abs(grid.x_grid(1) - grid.x_grid(2));
hv = abs(grid.v_grid(1) - grid.v_grid(2));
dx = abs(grid.poisson_grid(1) - grid.poisson_grid(2));
grid.domain_specs = [0, L, dx, 2; 0, L, hx, 2; -vmax, vmax, hv, 3];

%Total number of initial particles in computational domain
num_initial_particles = length(grid.x_grid)*length(grid.v_grid);

%Column 1 is Charge
%Column 2 is x positions
%Column 3 is v position
particle = NaN(num_initial_particles,5);
particle(:,2) = repmat(grid.x_grid, length(grid.v_grid), 1);
particle(:,3) = repelem(grid.v_grid,length(grid.x_grid));
particle(:,4) = ones(num_initial_particles, 1);
particle(:,5) = zeros(num_initial_particles,1);

switch problem
    case 'Linear_Landau_Damping' 
        
        initial_cond = @(x,v,alpha,k)(1/sqrt(2*pi)*exp(-v^2/2)*(1+alpha*cos(k*x)));
        
    case 'Two_Stream_Instability'
        
        initial_cond = @(x,v,alpha,k)(1/sqrt(2*pi)*v^2*exp(-v^2/2)*(1+alpha*cos(k*x)));
        
end

        
        

%Samples the distribution for IC's
for m = 1:num_initial_particles
    
    particle(m,1) = hx*hv*initial_cond(particle(m,2), particle(m,3), alpha, k);
 
end


%Removes all particles that don't have a weight of at least 10^-9
particle(particle(:,1) < 10^(-9),:) = [];


plotting_output = particle_push_MEX(particle, grid, T, dt, Remapping);


end