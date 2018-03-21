% Complexity Plots 

% Serial
nserial = [10000, 20000, 40000, 80000, 100000];
sim_time_serial = [1.39049, 2.92219, 5.99823, 12.7745, 28.5469];
slopeestimatesSerial = [1.071457, 1.037487, 1.090658, 1.160067];
slopeavgSerial = 1.084748;

% OpenMP 
nomp = [50000, 100000, 200000, 300000, 600000, 900000, 1200000, 1600000];
sim_time_omp = [3.17237, 6.62514, 14.6321, 27.5341, 78.303, 123.969, 171.594, 235.387];

polyfit(log(nomp), log(sim_time_omp), 1)

for j=2:length(nomp)
    slopeomp(j-1) = (log(sim_time_omp(j)) - log(sim_time_omp(j-1)))/(log(nomp(j)) - log(nomp(j-1)));
end

avgslopeomp = mean(slopeomp)


% loglog(nserial, sim_time_serial, 'b-o')
% set(gcf, 'color', 'w')
% xlabel('log(N)', 'Interpreter', 'LaTeX')
% ylabel('log(Time) (s)', 'Interpreter', 'LaTeX')
% title('Complexity', 'Interpreter', 'LaTeX')
% hold on
% loglog(nomp, sim_time_omp, 'r-o')
% legend('Serial', 'OpenMP')


% MPI 
nmpi = [100000, 200000, 300000, 600000, 900000, 1200000, 1600000];
sim_time_mpi = [0.61862, 1.25712, 1.94946, 4.00544, 5.39331, 8.33291, 11.344];

polyfit(log(nmpi), log(sim_time_mpi), 1)

for j=2:length(nmpi)
    slopempi(j-1) = (log(sim_time_mpi(j)) - log(sim_time_mpi(j-1)))/(log(nmpi(j)) - log(nmpi(j-1)));
end

avgslopempi = mean(slopempi)

% loglog(nmpi, sim_time_mpi, 'r-o')
% set(gcf, 'color', 'w')
% xlabel('log(N)', 'Interpreter', 'LaTeX')
% ylabel('log(Time) (s)', 'Interpreter', 'LaTeX')
% title('Complexity', 'Interpreter', 'LaTeX')
% hold on



particles = [50000      100000      200000      400000      800000]; 
sim_time_gpu = [0.407502, 0.9216, 1.96107, 4.08213, 8.35813];

for j=2:length(sim_time_gpu)
    slopegpu(j-1) = (log(sim_time_gpu(j)) - log(sim_time_gpu(j-1)))/(log(particles(j)) - log(particles(j-1)));
end

avgslopegpu = mean(slopegpu)

loglog(particles, sim_time_gpu, 'r-o')
set(gcf, 'color', 'w')
xlabel('log(# particles)', 'Interpreter', 'LaTeX')
ylabel('log(Time) (s)', 'Interpreter', 'LaTeX')
title('Complexity', 'Interpreter', 'LaTeX')
hold on
