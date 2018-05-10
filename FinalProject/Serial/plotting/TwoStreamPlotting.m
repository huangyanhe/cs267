%type LLD1p1pDTestEnergy
clear;
fileID = fopen('TwoStreamT15', 'r');
j=1;
%endline = 3774;
endline = 4014;
for i=1:endline
    c= fgetl(fileID);
    if i>=241
        %if i>=6
        val = str2num(c);
        if mod(i-241,3)==0
            x(j,1) = val;
        elseif mod(i-241,3)==1
            v(j,1) = val;
        else
            charge(j,1) = val;  
            j = j+1;
        end
        
    end
end
fclose(fileID);

[xi,yi] = meshgrid(0:0.1:4*pi, -10:0.1:10);
zi = griddata(x, v, charge,xi,yi);
s1 = surf(xi,yi,zi);
%scatter3(particle_out(j).values(:,2), particle_out(j).values(:,3), particle_out(j).values(:,1), 'bo')
axis([0 4*pi -10 10 -1 1])
view(0,90)
xlabel('position', 'Interpreter', 'LaTex')
ylabel('velocity', 'Interpreter', 'LaTex')
colorbar
s1.EdgeColor = 'none';
drawnow;
