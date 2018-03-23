function plotting_output = particle_push_MEX(particle, grid, T, dt, Remapping)


plotting_output(1).particle_time_series = particle;
plotting_output(1).electric_field_amplitude = [0, 0];
plotting_output(1).remap_times = NaN;
hx = grid.domain_specs(2,3);
dx = grid.domain_specs(1,3);
num_time_steps = floor(T/dt);

%%Main Time Loop
for it=1:num_time_steps
    
    if mod(it,5)==0 && Remapping == 1
    %if  max(abs(particle(:,4) - 1))>= 3/100*hx && Remapping == 1 
        (it-1)*dt
        plotting_output(it+1).remap_times = (it-1)*dt;
        particle = remapping_2d_MEX(particle, grid);
        particle(particle(:,1) < 10^(-9),:) = [];
        
    else
        
                plotting_output(it+1).remap_times = NaN;

        
    end
    
    [k1, k2, k3, particle, electric_field_amplitude, DE] = k_generator_MEX(particle, grid, dt);
    
    particle(:,2) = particle(:,2) + particle(:,3)*dt + 1/6*dt^2*(k1 + 2*k2);
    particle(:,3) = particle(:,3) + 1/6*dt*(k1 + 4*k2 + k3);
    particle(:,4) = particle(:,4) + dt*particle(:,5);
    particle(:,5) = particle(:,5) + dt*particle(:,4).*DE(:,1);
    
    particle = particle_wrapping(particle, grid);
    
    plotting_output(it+1).particle_time_series = particle;
    plotting_output(it+1).electric_field_amplitude = [(it-1)*dt, electric_field_amplitude];

end

%Final Time Step
dt = T-num_time_steps*dt;

[k1, k2, k3, particle, electric_field_amplitude, DE] = k_generator_MEX(particle, grid, dt);

particle(:,2) = particle(:,2) + particle(:,3)*dt + 1/6*dt^2*(k1 + 2*k2);
particle(:,3) = particle(:,3) + 1/6*dt*(k1 + 4*k2 + k3);

particle = particle_wrapping(particle, grid);

plotting_output(num_time_steps+2).particle_time_series = particle;
plotting_output(num_time_steps+2).electric_field_amplitude = [T, electric_field_amplitude];
plotting_output(num_time_steps+2).remap_times = NaN;

end