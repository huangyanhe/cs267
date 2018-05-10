% # threads = 4
% One node Scaling
% running 100 time steps 
n = [1       2       4      8]; 
M = 10;
times1 = [196.853 107.794 56.1449 30.238];
%serialTime = 642.68;

%Final time is 0.25. Time step changes so need to scale these
% Should double check this formula.
strongEfficiency1 = times1(1)./(n.*times1);

figure(1)
plot(n, strongEfficiency1, 'k-x')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Efficiency', 'Interpreter', 'LaTeX')
title('Strong Scaling', 'Interpreter', 'LaTeX')
axis([0 inf 0 1])
hold on 

% Multinode scaling
% running 100 time steps 
N = [2       4      8       16      32]; 
M = 10;
times2 = [15.0603    7.59036     4.43485     3.2736      2.9423];
%serialTime = 642.68;

%Final time is 0.25. Time step changes so need to scale these
% Should double check this formula.
strongEfficiency2 = times2(1)*2./(N.*times2);

plot(N, strongEfficiency2, 'r-x')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Efficiency', 'Interpreter', 'LaTeX')
title('Strong Scaling', 'Interpreter', 'LaTeX')
axis([0 inf 0 1])
hold on 
legend('Scaling MPI Processes Per Node', 'Scaling Number of Nodes')

%Number of MPI tasks
MPITasks1 = 2.^[0:1:3];
MPITasks2 = 2.^[4:1:8];
figure(2)

plot(MPITasks1, times1, 'k-x')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Efficiency', 'Interpreter', 'LaTeX')
title('Strong Scaling', 'Interpreter', 'LaTeX')
axis([0 inf 0 1])
hold on

plot(MPITasks2, times2, 'r-x')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('time', 'Interpreter', 'LaTeX')
title('Strong Scaling', 'Interpreter', 'LaTeX')
axis([0 inf 0 inf])
legend('Scaling MPI Processes Per Node', 'Scaling Number of Nodes')

figure(3)
loglog([MPITasks1 MPITasks2], [times1 times2], 'b-x')
set(gcf, 'color', 'w')
xlabel('$\log$(MPI Tasks)', 'Interpreter', 'LaTeX')
ylabel('$\log(t)$', 'Interpreter', 'LaTeX')
title('Strong Scaling LogLog plot', 'Interpreter', 'LaTeX')
axis([0 inf 0 inf])

hold on
loglog([MPITasks1 MPITasks2], times1(1)./[MPITasks1 MPITasks2], 'k--')

MPITasks = [MPITasks1 MPITasks2];
times = [times1 times2];
eff = times(1)./(MPITasks.* times);

figure(4)
plot(MPITasks, eff, 'k-x')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Efficiency', 'Interpreter', 'LaTeX')
title('Strong Scaling', 'Interpreter', 'LaTeX')
axis([0 inf 0 1])
hold on


