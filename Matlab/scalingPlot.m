% Efficiency Plots

% InitialCapacity = 2
processors = [1       2       4       6      12      18      24      32]; 
strongefficiancy = [0.93    0.88    0.63    0.60    0.56    0.51    0.51    0.47];
strongspeedup = [0.93    1.77    2.53    3.58    6.77    9.24   12.19   15.01];
weakefficiancy = [0.93    0.81    0.53    0.41    0.28    0.26    0.23    0.19];

strongspeedup = strongspeedup/strongspeedup(1);
% InitialCapacity = 5
%processors = [1       2       4       6      12      18      24      32]; 
%strongefficiancy = [0.94    0.87   0.58    0.57    0.51    0.44    0.41    0.44];
%strongspeedup = [0.94    1.74    2.34    3.43    6.16    7.98   9.78   14.12];
%weakefficiancy = [0.94    0.80    0.51   0.38    0.27    0.24    0.22    0.18];


figure(1)
subplot(1,3,1)
plot(processors, strongefficiancy, 'b-o')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Efficiency', 'Interpreter', 'LaTeX')
title('Strong Scaling', 'Interpreter', 'LaTeX')
axis([0 32 0 1])

%figure(2)
subplot(1,3,2)
plot(processors, strongspeedup, 'b-o')
hold on
plot(processors, processors, 'r-')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Speedup', 'Interpreter', 'LaTeX')
title('Strong Scaling Speedup', 'Interpreter', 'LaTeX')
axis([0 32 0 32])

%figure(3)
subplot(1,3,3)
plot(processors, weakefficiancy, 'b-o')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Efficiency', 'Interpreter', 'LaTeX')
title('Weak Scaling Efficiency', 'Interpreter', 'LaTeX')
axis([0 32 0 1])

clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MPI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
processors = [1       2       4       6      12      16      24      32]; 
strongefficiancy = [0.90    0.85    0.86    0.85    0.84    0.83    0.84    0.80];
strongspeedup = [0.90    1.70    3.43    5.09    10.04    13.35   20.13   25.47];
weakefficiancy = [0.90    0.82    0.78    0.76    0.72    0.72    0.71    0.70];

%comment out if you want the comparison to serial
% strongefficiancy = strongefficiancy./strongefficiancy(1);
% strongspeedup = strongspeedup./strongspeedup(1);
% weakefficiancy = weakefficiancy./weakefficiancy(1);
% sum(strongefficiancy)/8
% sum(weakefficiancy)/8



figure(2)
subplot(1,3,1)
plot(processors, strongefficiancy, 'b-o')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Efficiency', 'Interpreter', 'LaTeX')
title('Strong Scaling', 'Interpreter', 'LaTeX')
axis([0 32 0 1])

%figure(2)
subplot(1,3,2)
plot(processors, strongspeedup, 'b-o')
hold on
plot(processors, processors, 'r-')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Speedup', 'Interpreter', 'LaTeX')
title('Speedup', 'Interpreter', 'LaTeX')
axis([0 32 0 32])

%figure(3)
subplot(1,3,3)
plot(processors, weakefficiancy, 'b-o')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Efficiency', 'Interpreter', 'LaTeX')
title('Weak Scaling ', 'Interpreter', 'LaTeX')
axis([0 32 0 1])


particles = [50000      100000      200000      400000      800000]; 
timeGPU = [0.459329     1.01139     2.14141     4.42958     9.02934];

timeMPI12 = [0.944574     1.96794     5.4934     17.4693     38.4297];
timeMPI16 = [0.678022     1.36453     3.19875     11.4526     27.2849];
timeMPI20 = [0.530882     1.07822     2.56274     9.35543     22.0332];
timeMPI24 = [0.445845     0.897959     2.05022     7.78116     18.9094];
timeMPI28 = [0.391225    0.788692     2.44771     6.6163     16.5146];
timeMPI32 = [0.33126        0.692382    1.53979     5.73343     14.1358];

timeOMP12 = [1.66604    3.40598     7.34263     18.8619     40.1771];
timeOMP16 = [0.987774    1.93273     4.20508     12.9093     30.5746];
timeOMP20 = [0.742382    1.51571     5.08022     12.9279     29.6726];
timeOMP24 = [0.646606    1.30369     2.86733     8.49744     20.0645];
timeOMP28 = [0.562863    1.1487     2.50834     7.34085     17.5876];

figure(3)
plot(particles, timeGPU,'-o')
hold on
plot(particles, timeMPI12,'-o')
plot(particles, timeMPI16,'-o')
plot(particles, timeMPI20,'-o')
plot(particles, timeMPI24,'-o')
plot(particles, timeMPI28,'-o')
plot(particles, timeMPI32,'-o')
plot(particles, timeOMP12,'-o')
plot(particles, timeOMP16,'-o')
plot(particles, timeOMP20,'-o')
plot(particles, timeOMP24,'-o')
plot(particles, timeOMP28,'-o')
set(gcf, 'color', 'w')
xlabel('Particles', 'Interpreter', 'LaTeX')
ylabel('time', 'Interpreter', 'LaTeX')
title('GPU, MPI, OpenMP, ', 'Interpreter', 'LaTeX')
axis([50000 800000 -inf inf])
legend('GPU', 'MPI-12', 'MPI-16', 'MPI-20', 'MPI-24', 'MPI-28', 'MPI-32', 'OMP-12', 'OMP-16', 'OMP-20', 'OMP-24', 'OMP-28')


figure(4)
semilogx(particles, timeGPU,'-o')
hold on
semilogx(particles, timeMPI12,'-o')
semilogx(particles, timeMPI16,'-o')
semilogx(particles, timeMPI20,'-o')
semilogx(particles, timeMPI24,'-o')
semilogx(particles, timeMPI28,'-o')
semilogx(particles, timeMPI32,'-o')
semilogx(particles, timeOMP12,'-o')
semilogx(particles, timeOMP16,'-o')
semilogx(particles, timeOMP20,'-o')
semilogx(particles, timeOMP24,'-o')
semilogx(particles, timeOMP28,'-o')
set(gcf, 'color', 'w')
xlabel('Particles', 'Interpreter', 'LaTeX')
ylabel('time', 'Interpreter', 'LaTeX')
title('GPU, MPI, OpenMP, ', 'Interpreter', 'LaTeX')
axis([50000 800000 -inf inf])
legend('GPU', 'MPI-12', 'MPI-16', 'MPI-20', 'MPI-24', 'MPI-28', 'MPI-32', 'OMP-12', 'OMP-16', 'OMP-20', 'OMP-24', 'OMP-28')


sim_time_gpu = [0.407502, 0.9216, 1.96107, 4.08213, 8.35813];
sim_time_serialgpu = [8.69673, 21.3359, 46.0087, 145.731, 328.637];

figure(5)
plot(particles, sim_time_gpu,'-o')
hold on
plot(particles, sim_time_serialgpu,'-o')
set(gcf, 'color', 'w')
xlabel('Particles', 'Interpreter', 'LaTeX')
ylabel('time', 'Interpreter', 'LaTeX')
title('GPU, MPI, OpenMP, ', 'Interpreter', 'LaTeX')
axis([50000 800000 -inf inf])
legend('GPU', 'Serial')


