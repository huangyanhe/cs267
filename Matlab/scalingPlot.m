% Efficiency Plots

% InitialCapacity = 2
processors = [1       2       4       6      12      18      24      32]; 
strongefficiancy = [0.93    0.88    0.63    0.60    0.56    0.51    0.51    0.47];
strongspeedup = [0.93    1.77    2.53    3.58    6.77    9.24   12.19   15.01];
weakefficiancy = [0.93    0.81    0.53    0.41    0.28    0.26    0.23    0.19];

strongspeedup = strongspeedup/strongspeedup(1)
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



