% Seems like this is wrong 
% 
% % For small test problems calculated to t = 0.25. renormalized for 1 step
% % work.
% processors = [1       4       16      64]; 
% M = [3 4 5 6];
% times = [0.667746 0.393189 0.735356 1.45];
% dt = 2./(2.^M);
% numTimeSteps = 0.25./dt; 
% serialTime = 0.405715;
% 
% %Final time is 0.25. Time step changes so need to scale these
% % Should double check this formula.
% %weakEfficiency = (serialTime/1/4)./(times./numTimeSteps);
% 
% 
% figure(1)
% plot(processors, weakEfficiency, 'b-o')
% hold on
% set(gcf, 'color', 'w')
% xlabel('Processors', 'Interpreter', 'LaTeX')
% ylabel('Efficiency', 'Interpreter', 'LaTeX')
% title('Weak Scaling', 'Interpreter', 'LaTeX')
% axis([0 64 0 1])


% Runs for larger test probems starting at M=10. Run for 100 steps.
processorsLarge = [1       4      16      64]; 
MLarge = [10 11 12];
timesLarge = [1042.26 962.985 865.35 764.034];
serialTimeLarge = 642.68;

%Final time is 0.25. Time step changes so need to scale these
% Should double check this formula.
weakEfficiency = (serialTimeLarge)./(timesLarge);


plot(processorsLarge, weakEfficiency, 'r-o')
set(gcf, 'color', 'w')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Efficiency', 'Interpreter', 'LaTeX')
title('Weak Scaling', 'Interpreter', 'LaTeX')
axis([0 inf 0 1])
hold on