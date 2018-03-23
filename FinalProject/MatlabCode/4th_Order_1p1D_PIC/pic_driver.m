%close all;
need_to_run_simulation = 1;

if need_to_run_simulation == 1
    %clear all;
    mex remapping.c
    mex deposition.c
    problem = 'Linear_Landau_Damping';
    Grid_Size = 8;
    Remapping = 1;
    tic
    plotting_output = pic_sub_driver(problem, Grid_Size, Remapping);
    toc
end
plotting = 0;


switch problem
    
    case 'Linear_Landau_Damping'
        
        size_data_out = size(plotting_output);
        electric_energy = NaN(size_data_out(:,2), 2);
        remap_markers = NaN(size_data_out(:,2), 2);
        
        for j=1:size_data_out(:,2)
            
            electric_energy(j, :) = plotting_output(j).electric_field_amplitude;
            remap_markers(j,:) = [plotting_output(j).remap_times, electric_energy(j, 2)];
            
        end
        
        figure(1);
        semilogy(electric_energy(:,1), 0.011*exp(-0.1533*electric_energy(:,1)), 'r--')
        hold on
        semilogy(electric_energy(:,1), electric_energy(:,2), 'k-');
        xlabel('time')
        ylabel('$\log(\xi_E)$', 'Interpreter', 'LaTex')
        title_str = sprintf('Grid Size = %f', Grid_Size);
        title(title_str, 'Interpreter', 'LaTex')
        set(gcf,'color','w')
        
        
        
        if plotting == 1
        
        make_movie = 0;
        [xi,yi] = meshgrid(0:0.1:4*pi, -10:0.1:10);

        figure(2);
        
        for j=1:size_data_out(:,2)
            
            %figure(1)
            subplot(1,3,1)
            zi = griddata(plotting_output(j).particle_time_series(:,2), plotting_output(j).particle_time_series(:,3), plotting_output(j).particle_time_series(:,1),xi,yi);
            s1 = surf(xi,yi,zi);
            %scatter3(particle_out(j).values(:,2), particle_out(j).values(:,3), particle_out(j).values(:,1), 'bo')
            axis([0 4*pi -10 10 -1 1])
            view(0,90)
            xlabel('position', 'Interpreter', 'LaTex')
            ylabel('velocity', 'Interpreter', 'LaTex')
            colorbar      
            s1.EdgeColor = 'none';
            drawnow;

            
            %figure(2)
            subplot(1,3,2)
            zi2 = griddata(plotting_output(j).particle_time_series(:,2), plotting_output(j).particle_time_series(:,3), plotting_output(j).particle_time_series(:,4),xi,yi);
            s2 = surf(xi,yi,zi2);
            %scatter3(plotting_output(j).particle_time_series(:,2), plotting_output(j).particle_time_series(:,3), plotting_output(j).particle_time_series(:,4), 'bo')
            axis([0 4*pi -10 10])
            view(0,90)
            xlabel('position', 'Interpreter', 'LaTex')
            colorbar
            set(gcf,'color','w')
            s2.EdgeColor = 'none';
            title('F', 'Interpreter', 'LaTex')
            drawnow;

                        
            %figure(3)
            subplot(1,3,3)
            zi3 = griddata(plotting_output(j).particle_time_series(:,2), plotting_output(j).particle_time_series(:,3), plotting_output(j).particle_time_series(:,5),xi,yi);
            s3 = surf(xi,yi,zi3);
            %scatter3(plotting_output(j).particle_time_series(:,2), plotting_output(j).particle_time_series(:,3), plotting_output(j).particle_time_series(:,5), 'bo')
            axis([0 4*pi -10 10])
            view(0,90)
            xlabel('position', 'Interpreter', 'LaTex')
            colorbar
            set(gcf,'color','w')
            s3.EdgeColor = 'none';
            title('G', 'Interpreter', 'LaTex')
            colorbar
            time = num2str(j*1/Grid_Size,3);
            B = uicontrol('Style', 'text',...
                'String', ['t=' time ],... %replace something with the text you want
                'Units','normalized',...
                'Position', [0.9 0.9 0.1 0.1]);
            set(B,'backgroundcolor', 'w')
            drawnow;
            
        end
        
        end
        
    case 'Two_Stream_Instability'

        if plotting == 1
        
        make_movie = 0;
        [xi,yi] = meshgrid(0:0.1:4*pi, -10:0.1:10);
        
        if make_movie == 1
            v = VideoWriter('Two_Stream_FG_Tracking.mp4','MPEG-4');
            open(v);
        end
        
        size_data_out = size(plotting_output);
        
        figure;
        
        for j=1:size_data_out(:,2)
            
            %figure(1)
            subplot(1,3,1)
            zi = griddata(plotting_output(j).particle_time_series(:,2), plotting_output(j).particle_time_series(:,3), plotting_output(j).particle_time_series(:,1),xi,yi);
            s1 = surf(xi,yi,zi);
            %scatter3(particle_out(j).values(:,2), particle_out(j).values(:,3), particle_out(j).values(:,1), 'bo')
            axis([0 4*pi -10 10 -1 1])
            view(0,90)
            xlabel('position', 'Interpreter', 'LaTex')
            ylabel('velocity', 'Interpreter', 'LaTex')
            colorbar      
            s1.EdgeColor = 'none';
            drawnow;

            
            %figure(2)
            subplot(1,3,2)
            zi2 = griddata(plotting_output(j).particle_time_series(:,2), plotting_output(j).particle_time_series(:,3), plotting_output(j).particle_time_series(:,4),xi,yi);
            s2 = surf(xi,yi,zi2);
            %scatter3(plotting_output(j).particle_time_series(:,2), plotting_output(j).particle_time_series(:,3), plotting_output(j).particle_time_series(:,4), 'bo')
            axis([0 4*pi -10 10])
            view(0,90)
            xlabel('position', 'Interpreter', 'LaTex')
            colorbar
            set(gcf,'color','w')
            s2.EdgeColor = 'none';
            title('F', 'Interpreter', 'LaTex')
            drawnow;

                        
            %figure(3)
            subplot(1,3,3)
            zi3 = griddata(plotting_output(j).particle_time_series(:,2), plotting_output(j).particle_time_series(:,3), plotting_output(j).particle_time_series(:,5),xi,yi);
            s3 = surf(xi,yi,zi3);
            %scatter3(plotting_output(j).particle_time_series(:,2), plotting_output(j).particle_time_series(:,3), plotting_output(j).particle_time_series(:,5), 'bo')
            axis([0 4*pi -10 10])
            view(0,90)
            xlabel('position', 'Interpreter', 'LaTex')
            colorbar
            set(gcf,'color','w')
            s3.EdgeColor = 'none';
            title('G', 'Interpreter', 'LaTex')
            colorbar
            time = num2str(j*1/Grid_Size,3);
            B = uicontrol('Style', 'text',...
                'String', ['t=' time ],... %replace something with the text you want
                'Units','normalized',...
                'Position', [0.9 0.9 0.1 0.1]);
            set(B,'backgroundcolor', 'w')
            drawnow;
            
            if make_movie == 1
                
                frame=getframe(gcf);
                writeVideo(v, frame);
                
            end
            
        end
        
        
        
        if make_movie == 1
            close(v);
        end
        
        end
        
end