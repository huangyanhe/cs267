%type LLD1p1pDTestEnergy
clear;
fileID = fopen('LLD1p1DTestEnergyM4W40', 'r');
j=0;
endline = 245;
%endline = 485;
%endline = 965;
for i=1:endline
    c= fgetl(fileID);
    if i>=6
        val = str2num(c);
        if mod(i,2)==1
            time(j,1) = val;
        else
            j = j+1;
            electric_energy(j,1) = val;
        end
        
    end
end
fclose(fileID);

%semilogy(time(:,1), 0.011*exp(-0.1533*time(:,1)), 'r--')
semilogy(time(:,1), electric_energy(1,1)*exp(-0.394*time(:,1)), 'r--')
hold on
semilogy(time(:,1), electric_energy(:,1), 'k-');
xlabel('time')
ylabel('$\log(\xi_E)$', 'Interpreter', 'LaTex')
%title_str = sprintf('Grid Size = %f', Grid_Size);
%title(title_str, 'Interpreter', 'LaTex')
set(gcf,'color','w')