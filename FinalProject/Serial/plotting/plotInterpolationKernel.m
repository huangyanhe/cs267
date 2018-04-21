load ../src/InterpolationKernel/UnitTests/plotInterpolationKernelW4r0


plot(plotInterpolationKernelW4r0(:,1), plotInterpolationKernelW4r0(:,2))
title('Interpolation Kernel', 'Interpreter', 'LaTex')
xlabel('x', 'Interpreter', 'LaTex')
ylabel('W(x)', 'Interpreter', 'LaTex')
set(gcf,'color','w')