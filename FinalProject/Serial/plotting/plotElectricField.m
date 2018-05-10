%type LLD1p1pDTestEnergy
clear;
fileID = fopen('LLD1p1DTestEnergyM6W20', 'r');
j=1;
%endline = 245;
%endline = 485;
%endline = 965;
%endline = 3840;
endline = 1920
for i=1:endline
    c= fgetl(fileID);
    if i>=1
    %if i>=6
        val = str2num(c);
        if mod(i,2)==0
            time(j,1) = val;
        else
            j = j+1;
            electric_energy(j,1) = val;
        end
        
    end
end
fclose(fileID);

semilogy(time(:,1), 0.011*exp(-0.1533*time(:,1)), 'r--')
%semilogy(time(:,1), electric_energy(2)*exp(-0.395*time(:,1)), 'r--')
hold on
semilogy(time(:,1), electric_energy(:,1), 'k-');
xlabel('time')
ylabel('$\log(\xi_E)$', 'Interpreter', 'LaTex')
%title_str = sprintf('Grid Size = %f', Grid_Size);
%title(title_str, 'Interpreter', 'LaTex')
set(gcf,'color','w')