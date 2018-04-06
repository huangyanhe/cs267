function phi = poisson_solver_FFT(density, grid)

%density = density - mean(density);
num_points_x = length(grid.poisson_grid);
L = grid.domain_specs(1,2);

Nx = num_points_x;

transDensity = fft(density);

k = [0:Nx/2 -Nx/2+1:-1];
%k=[0:Nx/2-1 0 -Nx/2+1:-1]
kx = ((2*pi/L.*k).^2)';
%kx
%transDensity(1) = 0;
transDensity(2:end) = transDensity(2:end)./kx(2:end);
%transDensity
%ifft(transDensity)

%phi = -real(ifft(transDensity));
phi = ifft(transDensity, 'symmetric');

end