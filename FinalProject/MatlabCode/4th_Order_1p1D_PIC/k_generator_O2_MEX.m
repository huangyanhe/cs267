function [k1, k2, particle, electric_field_amplitude, kF1, kF2, kG1, kG2] = k_generator_O2_MEX(particle, grid, dt)
    
    dx = grid.domain_specs(1,3);
    L = grid.domain_specs(1,2);
    
    %% Compute k1
    density = density_remapping_MEX(particle, grid);
    phi = poisson_solver(density, grid);
    e_field = gradient(phi, grid);
    forces = force_remapping_MEX(particle, grid, e_field);
    k1 = forces;
    
    electric_field_amplitude = sqrt(dx/L*(transpose(e_field)*e_field));
    
    %% Compute k2
    particle2 = particle;
    particle2(:,2) = particle2(:,2) + dt*particle2(:,3);
    particle2 = particle_wrapping(particle2, grid);
    
    density2 = density_remapping_MEX(particle2, grid);
    phi2 = poisson_solver(density2, grid);
    e_field2 = gradient(phi2, grid);
    forces2 = force_remapping_MEX(particle2, grid, e_field2);
    k2 = forces2;
        
    
    %% Compute kFs and kGs
    dE1 = -gradient(e_field, grid);
    dE1_k = force_remapping_MEX(particle, grid, dE1);
    
    particle= particle2;
   
    dE2 = -gradient(e_field2, grid);
    dE2_k = force_remapping_MEX(particle, grid, dE2);
    
    kF1 = particle(:, 5);
    kG1 = dE1_k.*particle(:,4);
    kF2 = particle(:,5) + dt*kG1;
    kG2 = dE2_k.*(particle(:, 4) + dt*kF1);
    
    

end