function particle = particle_wrapping(particle, grid)
    
    L = grid.domain_specs(1,2); 
    TE = grid.domain_specs(3,2);
    BE = grid.domain_specs(3,1);
    Lv = abs(TE - BE);
    
    %Imposes the periodicity of the boundary conditions
    particle(particle(:,2) >= L, 2) = particle(particle(:,2) >= L, 2) - L;
    particle(particle(:,2) < 0, 2) = particle(particle(:,2) < 0, 2) + L;

    %Imposes the periodicity of the boundary conditions
    particle(particle(:,3) >= TE, 3) = particle(particle(:,3) >= TE, 3) - Lv;
    particle(particle(:,3) < BE, 3) = particle(particle(:,3) < BE, 3) + Lv;

end