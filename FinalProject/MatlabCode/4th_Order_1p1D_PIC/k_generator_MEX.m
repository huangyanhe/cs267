function [k1, k2, k3, particle, electric_field_amplitude, DE] = k_generator_MEX(particle, grid, dt)
    
    dx = grid.domain_specs(1,3);
    L = grid.domain_specs(1,2);
    
    %% Compute k1
    density = density_remapping_MEX(particle, grid);
    phi = poisson_solver_FFT(density, grid);
    e_field = gradient(phi, grid);
    forces = force_remapping_MEX(particle, grid, e_field);
    k1 = forces;
    
    electric_field_amplitude = sqrt(dx/L*(transpose(e_field)*e_field));
    
    %% Compute k2
    particle2 = particle;
    particle2(:,2) = particle2(:,2) + 1/2*dt*particle2(:,3) +1/8*k1*dt^2;
    particle2 = particle_wrapping(particle2, grid);
    
    density2 = density_remapping_MEX(particle2, grid);
    phi2 = poisson_solver_FFT(density2, grid);
    e_field2 = gradient(phi2, grid);
    forces2 = force_remapping_MEX(particle2, grid, e_field2);
    k2 = forces2;
    
    %% Compute k3
    particle3 = particle;
    particle3(:,2) = particle3(:,2) + dt*particle3(:,3) +1/2*k2*dt^2;
    particle3 = particle_wrapping(particle3, grid);
    
    density3 = density_remapping_MEX(particle3, grid);
    phi3 = poisson_solver_FFT(density3, grid);
    e_field3 = gradient(phi3, grid);
    forces3 = force_remapping_MEX(particle3, grid, e_field3);
    k3 = forces3;
    
    %% Compute kFs and kGs
    dE1 = -gradient(e_field, grid);
    DE(:,1) = force_remapping_MEX(particle, grid, dE1);
        
    
    

end