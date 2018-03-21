nodes = [1      2      4      8]; 
timeHumanTotal = [22.307675     31.252265     23.411704     13.792652];
timeHumanInsert = [16.317191     20.029628     14.587123     7.844210];
timeLargeTotal = [6.918531     9.540548     7.101898     4.006188];
timeLargeInsert = [4.993300     5.572123     3.817618    2.340963];
timeTestTotal = [1.116526     1.438917     1.309812     0.776025];
timeTestInsert = [0.789641      0.911556     0.834655     0.482836];

figure(1)
plot(nodes, timeHumanTotal,'-o')
hold on
plot(nodes, timeLargeTotal,'-o')
plot(nodes, timeTestTotal,'-o')
plot(nodes, timeHumanInsert,'-o')
plot(nodes, timeLargeInsert,'-o')
plot(nodes, timeTestInsert,'-o')
xlabel('Nodes (w/32 Processors/Node)', 'Interpreter', 'LaTeX')
ylabel('time', 'Interpreter', 'LaTeX')
title('Node Scaling', 'Interpreter', 'LaTeX')
axis([1 8 0 35])
legend('Human Total', 'Large Total', 'Test Total', 'Human Insertion', 'Large Insertion', 'Test Insertion')
set(gcf, 'color', 'w')

Processors = [1     2      4      8   16    32    64]; 
timeHumanTotalP = [235.729162      199.620335     124.953682     70.877791     41.494250     22.670032       19.971469];
timeHumanInsertP = [119.374428      137.511313     91.382668     52.296240     29.624552     16.628471       13.383914];
timeLargeTotalP = [72.483718     59.794885     36.041603    21.277969      15.673672        8.670716     7.514075];
timeLargeInsertP = [35.782025    41.186960     26.087343    15.621722      11.534327        6.228951     4.828977];
timeTestTotalP = [17.292061     9.436486     8.242512     3.511923      2.635697        1.113242     1.036404];
timeTestInsertP = [8.902240    6.563542     5.858926    2.615434      1.554115         0.756301    0.658039];

figure(2)
plot(Processors, timeHumanTotalP,'-o')
hold on
plot(Processors, timeLargeTotalP,'-o')
plot(Processors, timeTestTotalP,'-o')
plot(Processors, timeHumanInsertP,'-o')
plot(Processors, timeLargeInsertP,'-o')
plot(Processors, timeTestInsertP,'-o')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('time', 'Interpreter', 'LaTeX')
title('Processor Scaling (N=1)', 'Interpreter', 'LaTeX')
axis([1 64 0 250])
legend('Human Total', 'Large Total', 'Test Total', 'Human Insertion', 'Large Insertion', 'Test Insertion')
set(gcf, 'color', 'w')



figure(3)
loglog(Processors, timeHumanTotalP,'-o')
hold on
loglog(Processors, timeLargeTotalP,'-o')
loglog(Processors, timeTestTotalP,'-o')
loglog(Processors, timeHumanInsertP,'-o')
loglog(Processors, timeLargeInsertP,'-o')
loglog(Processors, timeTestInsertP,'-o')
xlabel('$\log$(Processors)', 'Interpreter', 'LaTeX')
ylabel('$\log$(time)', 'Interpreter', 'LaTeX')
title('Processor Scaling (N=1)', 'Interpreter', 'LaTeX')
axis([1 64 0 250])
legend('Human Total', 'Large Total', 'Test Total', 'Human Insertion', 'Large Insertion', 'Test Insertion')
set(gcf, 'color', 'w')