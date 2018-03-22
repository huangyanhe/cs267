nodes = [1      2      4      8]; 
timeHumanTotal = [22.307675     31.252265     23.411704     13.792652];
timeHumanInsert = [16.317191     20.029628     14.587123     7.844210];
timeHumanRead = [5.768747     10.344762     6.752962     3.098957];

timeLargeTotal = [6.918531     9.540548     7.101898     4.006188];
timeLargeInsert = [4.993300     5.572123     3.817618    2.340963];
timeLargeRead = [1.658085     2.169158     1.532800    0.770053];

timeTestTotal = [1.116526     1.438917     1.309812     0.776025];
timeTestInsert = [0.789641      0.911556     0.834655     0.482836];
timeTestRead = [0.245904      0.368123     0.373610      0.148624];

figure(1)
plot(nodes, timeHumanTotal,'-o')
hold on
plot(nodes, timeLargeTotal,'-o')
plot(nodes, timeTestTotal,'-o')
plot(nodes, timeHumanInsert,'-o')
plot(nodes, timeLargeInsert,'-o')
plot(nodes, timeTestInsert,'-o')
plot(nodes, timeHumanRead,'-o')
plot(nodes, timeLargeRead,'-o')
plot(nodes, timeTestRead,'-o')
xlabel('Nodes (w/32 Processors/Node)', 'Interpreter', 'LaTeX')
ylabel('time', 'Interpreter', 'LaTeX')
title('Node Scaling', 'Interpreter', 'LaTeX')
axis([1 8 0 35])
legend('Human Total', 'Large Total', 'Test Total', 'Human Insertion', 'Large Insertion', 'Test Insertion', 'Human Read', 'Large Read', 'Test Read')
set(gcf, 'color', 'w')

Processors = [1     2      4      8   16    32    64]; 
timeHumanTotalP = [235.729162      199.620335     124.953682     70.877791     41.494250     22.670032       19.971469];
timeHumanInsertP = [119.374428      137.511313     91.382668     52.296240     29.624552     16.628471       13.383914];
timeHumanReadP = [116.354146      62.108153     33.570244     17.642998     10.421114     5.874709       5.995908];
timeLargeTotalP = [72.483718     59.794885     36.041603    21.277969      15.673672        8.670716     7.514075];
timeLargeInsertP = [35.782025    41.186960     26.087343    15.621722      11.534327        6.228951     4.828977];
timeLargeReadP = [36.700924    18.607071     9.727755    5.237618      3.713494      2.033984     2.158916];
timeTestTotalP = [17.292061     9.436486     8.242512     3.511923      2.635697        1.113242     1.036404];
timeTestInsertP = [8.902240    6.563542     5.858926    2.615434      1.554115         0.756301    0.658039];
timeTestReadP = [8.389075    2.661325    1.988651    0.768317      0.504736        0.245379    0.292481];

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
loglog(Processors, timeHumanReadP,'-o')
loglog(Processors, timeLargeReadP,'-o')
loglog(Processors, timeTestReadP,'-o')
xlabel('$\log$(Processors)', 'Interpreter', 'LaTeX')
ylabel('$\log$(time)', 'Interpreter', 'LaTeX')
title('Processor Scaling (N=1)', 'Interpreter', 'LaTeX')
axis([1 64 0 250])
legend('Human Total', 'Large Total', 'Test Total', 'Human Insertion', 'Large Insertion', 'Test Insertion', 'Human Read', 'Large Read', 'Test Read')
set(gcf, 'color', 'w')


figure(4)
plot(Processors, 1./(Processors.*timeHumanTotalP/timeHumanTotalP(1)),'-o')
hold on
plot(Processors, 1./(Processors.*timeLargeTotalP/timeLargeTotalP(1)),'-o')
plot(Processors, 1./(Processors.*timeTestTotalP/timeTestTotalP(1)),'-o')
plot(Processors, 1./(Processors.*timeHumanInsertP/timeHumanInsertP(1)),'-o')
plot(Processors, 1./(Processors.*timeLargeInsertP/timeLargeInsertP(1)),'-o')
plot(Processors, 1./(Processors.*timeTestInsertP/timeTestInsertP),'-o')
plot(Processors, 1./(Processors.*timeHumanReadP/timeHumanReadP(1)),'-o')
plot(Processors, 1./(Processors.*timeLargeReadP/timeLargeReadP(1)),'-o')
plot(Processors, 1./(Processors.*timeTestReadP/timeTestReadP(1)),'-o')
xlabel('Processors', 'Interpreter', 'LaTeX')
ylabel('Efficiency', 'Interpreter', 'LaTeX')
title('Processor Scaling Efficiency (N=1)', 'Interpreter', 'LaTeX')
%axis([1 64 0 250])
legend('Human Total', 'Large Total', 'Test Total', 'Human Insertion', 'Large Insertion', 'Test Insertion', 'Human Read', 'Large Read', 'Test Read')
set(gcf, 'color', 'w')

figure(5)
plot(nodes, 1./(nodes.*timeHumanTotal/timeHumanTotal(1)),'-o')
hold on
plot(nodes, 1./(nodes.*timeLargeTotal/timeLargeTotal(1)),'-o')
plot(nodes, 1./(nodes.*timeTestTotal/timeTestTotal(1)),'-o')
plot(nodes, 1./(nodes.*timeHumanInsert/timeHumanInsert(1)),'-o')
plot(nodes, 1./(nodes.*timeLargeInsert/timeLargeInsert(1)),'-o')
plot(nodes, 1./(nodes.*timeTestInsert/timeTestInsert(1)),'-o')
plot(nodes, 1./(nodes.*timeHumanRead/timeHumanRead(1)),'-o')
plot(nodes, 1./(nodes.*timeLargeRead/timeLargeRead(1)),'-o')
plot(nodes, 1./(nodes.*timeTestRead/timeTestRead(1)),'-o')
xlabel('Nodes (w/32 Processors/Node)', 'Interpreter', 'LaTeX')
ylabel('Efficiency', 'Interpreter', 'LaTeX')
title('Node Scaling Efficiency', 'Interpreter', 'LaTeX')
%axis([1 8 0 35])
legend('Human Total', 'Large Total', 'Test Total', 'Human Insertion', 'Large Insertion', 'Test Insertion', 'Human Read', 'Large Read', 'Test Read')
set(gcf, 'color', 'w')
