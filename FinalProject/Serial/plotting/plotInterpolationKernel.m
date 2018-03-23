load ../src/InterpolationKernel/UnitTests/plotInterpolationKernelW0r0


plot(plotInterpolationKernelW0r0(:,1), plotInterpolationKernelW0r0(:,2))
title('Interpolation Kernel', 'Interpreter', 'LaTex')
xlabel('x', 'Interpreter', 'LaTex')
ylabel('W(x)', 'Interpreter', 'LaTex')
set(gcf,'color','w')