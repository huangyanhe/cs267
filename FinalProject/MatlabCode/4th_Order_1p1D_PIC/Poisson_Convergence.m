clear all; 
Cells = [4, 8, 16, 32, 64, 128];

errors = zeros(1, length(Cells));
h_vector = zeros(1, length(Cells));

for j = 1:length(Cells)
vmax =0; 
L = 1;

NumCells = Cells(j);
Numx = 8;
Numv = 0;

%Spatial Grid
spatial_grid_init = transpose(linspace(0,L,NumCells+1));
grid.poisson_grid = spatial_grid_init(1:end-1, :);


particle_spatial_grid_init = transpose(linspace(0,L,Numx+1));
grid.x_grid = particle_spatial_grid_init(1:end-1, :);
grid.v_grid = 0;


hx = abs(grid.x_grid(1) - grid.x_grid(2));
hv = 0;
dx = abs(grid.poisson_grid(1) - grid.poisson_grid(2));
h_vector(j) = dx;
grid.domain_specs = [0, L, dx, NaN; 0, L, hx, 2; -vmax, vmax, hv, 3];

density = (2*pi/L)^2*sin(2*pi*grid.poisson_grid/L);

phi = poisson_solver_FFT(density, grid);

plot(grid.poisson_grid, phi, 'rx')
hold on
plot(grid.poisson_grid, sin(2*pi*grid.poisson_grid/L), 'bo')
clf


errors(j) = max(abs(phi - sin(2*pi*grid.poisson_grid/L)));

end

loglog(h_vector, errors)
set(gcf, 'color', 'w')
xlabel('log(h)')
ylabel('log(errors)')
fit = polyfit(log(h_vector), log(errors), 1);
slope = fit(:,1)


