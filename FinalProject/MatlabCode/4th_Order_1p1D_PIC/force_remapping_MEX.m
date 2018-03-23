function forces = force_remapping_MEX(particle, grid, forces)
    
    %% Initialization for the For loops
    projected_particle = zeros(length(grid.poisson_grid), 2);
    projected_particle(:,2) = grid.poisson_grid;
    projected_particle(:,1) = forces;
    
    
    size_of_extension = 2;
    
    dx = grid.domain_specs(1,3);
    h_vector = dx;
    
    projected_particle = periodic_wrapping(projected_particle, h_vector, grid.domain_specs(1,:), size_of_extension);
    
    
    forces = deposition(projected_particle(:,1)', 1/dx.*projected_particle(:,2)', 1/dx.*particle(:,2)')';
    
    
    
    


end