function density = density_remapping_MEX(particle, grid)
    
    %% Initialization for the For loops
    
    
    size_of_extension = 2;
    
    h = grid.domain_specs(1,3);
    h_vector = h;
    
    particle = periodic_wrapping(particle, h_vector, grid.domain_specs(2,:), size_of_extension);
    
    
    density = deposition(1/h.*particle(:,1)', 1/h.*particle(:,2)', 1/h.*grid.poisson_grid')';

    


end