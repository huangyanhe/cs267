function remapped_particle = remapping_2d_MEX(particle_one_domain, grid)


%% Unwrapping of grid data structure
x_grid = grid.x_grid;
v_grid = grid.v_grid;
domain_specs = grid.domain_specs;


hx = domain_specs(2, 3);
hv = domain_specs(3 ,3);

num_remapped_particles = length(v_grid)*length(x_grid);

%% Initialization for the For loops
remapped_particle = zeros(num_remapped_particles, 5);
remapped_particle(:,2)= repelem(x_grid, length(v_grid));
remapped_particle(:,3) = repmat(v_grid, length(x_grid),1);
remapped_particle(:,4) = ones(num_remapped_particles, 1);
remapped_particle(:,5) = zeros(num_remapped_particles, 1);

%% Added stuff for dealing with periodic boundary conditions
size_of_extension = 3;
h_vector = [hx, hv];
particle = periodic_wrapping(particle_one_domain, h_vector, domain_specs(2:3,:), size_of_extension);


%% C++ Function Call
% Input: particle, remapped_particle, num_particles, num_remapped_particles
% Output: remapped_particle
remapped_particle(:,1) = remapping( particle(:,1)', 1/hx*particle(:,2)', 1/hv*particle(:,3)', 1/hx*remapped_particle(:,2)', 1/hv*remapped_particle(:,3)');


%% Positivity Preservation
remapped_particle = positivity_preservation(remapped_particle, grid);


end