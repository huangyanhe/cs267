function phi = poisson_solver(density, grid)

num_points_x = length(grid.poisson_grid);
delta_x = grid.domain_specs(1,3);

%Creates Spares 2nd order FD scheme for phi_0 to phi_(n-1) since we have
%perioidic boundary conditions
A = 1/(12*delta_x^2).*gallery('circul', [-30,16, -1,zeros([1, num_points_x-5]), -1, 16]);


%Enforces the periodic boundary conditions and makes the matrix nonsingular
%by requiring that the average is 0 
%I am not totally sure that this is the correct thing to do but the idea is
%taken from a comp. plasma phys. notes
A = [A; ones(1,length(A))*delta_x];
%density(end) = 0;
density = [density; 0];

phi = -A\density;


end