function e_field = gradient(phi, grid)

num_points_x = length(grid.poisson_grid);
delta_x = grid.domain_specs(1,3);

%Creates a sparse 1st derivative centered difference scheme for computing
%the electric field from the potential
B = 1/(12*delta_x).*gallery('circul', [0,8, -1,zeros([1, num_points_x-5]), 1, -8]);


e_field = -B*phi;

end