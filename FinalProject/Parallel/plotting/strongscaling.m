% MPI Scaling
% running 100 time steps 
processors = [1       4       16      64 256]; 
M = 10;
times = [1041.46 260.864 63.4973 17.4364 13.5267];
serialTime = 642.68;

%Final time is 0.25. Time step changes so need to scale these
% Should double check this formula.
strongEfficiency = serialTime./(processors.*times);

plot(processors, strongEfficiency, 'k-x')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Efficiency', 'Interpreter', 'LaTeX')
title('Strong Scaling', 'Interpreter', 'LaTeX')
axis([0 inf 0 1])
hold on 


% OMP Scaling
serialTime = 642.68;
processors = [1    2   4    8   16     32]; 
times = [650.765 342.914 197.027 126.717 89.1005 78.6867 ];
strongEfficiency = serialTime./(processors.*times);



plot(processors, strongEfficiency, 'r-x')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Efficiency', 'Interpreter', 'LaTeX')
title('Scaling', 'Interpreter', 'LaTeX')
axis([0 inf 0 1])