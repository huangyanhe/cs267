load ../src/InterpolationKernel/UnitTests/plotInterpolationKernelW2r0


plot(plotInterpolationKernelW2r0(:,1), plotInterpolationKernelW2r0(:,2))
title('Interpolation Kernel', 'Interpreter', 'LaTex')
xlabel('x', 'Interpreter', 'LaTex')
ylabel('W(x)', 'Interpreter', 'LaTex')
set(gcf,'color','w')