% running 100 time steps 
processors = [1       4       16      ]; 
M = 10;
times = [1041.46 260.864 63.4973 ];
serialTime = 642.68;

%Final time is 0.25. Time step changes so need to scale these
% Should double check this formula.
strongEfficiency = serialTime./(processors.*times);


figure(1)
plot(processors, strongEfficiency, 'b-o')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Efficiency', 'Interpreter', 'LaTeX')
title('Strong Scaling', 'Interpreter', 'LaTeX')
axis([0 inf 0 1])