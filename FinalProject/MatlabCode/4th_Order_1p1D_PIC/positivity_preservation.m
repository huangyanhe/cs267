function particle = positivity_preservation(particle, grid)


    hx = grid.domain_specs(2, 3);
    hv = grid.domain_specs(3 ,3);
    
    n = 0;

   while any(particle(:,1)< -10^(-14))
      
      negative_cell_index = find(particle(:,1)<0);
            
      particles_positive = particle(:,1);
      
       for j = 1:length(negative_cell_index) 
           
           index = negative_cell_index(j);
           
           index_neighbors = find(abs(particle(:,2) - particle(index,2))-hx<= 10^(-12)...
               & abs(particle(:,2) - particle(index,2))>=0 ...
               &  abs(particle(:,3) - particle(index,3))-hv<= 10^(-12)...
               &  abs(particle(:,3) - particle(index,3))>=0);
           
           %The point is a boundary point
           if length(index_neighbors)<9
               
               TE = grid.domain_specs(3,2);
               BE = grid.domain_specs(3,1);
               Lv = abs(TE - BE);
               LE = grid.domain_specs(2,1);
               RE = grid.domain_specs(2,2);
               Lx = abs(LE - RE);
               
               %Point is on an edge in v-space
               if abs(particle(index,3) - (TE - hv))<= 10^(-12) || abs(particle(index,3) - BE)<=10^(-12)
                   
                   index_wrapped1 = find(abs(particle(:,2) - particle(index,2))-hx<= 10^(-12)...
                       & abs(particle(:,2) - particle(index,2))>=0 ...
                       &  abs(particle(:,3) - particle(index,3)) - Lv<= 10^(-12)...
                       &  abs(particle(:,3) - particle(index,3))-(Lv-hv)>=-10^(-12) );
                   
               else
                   index_wrapped1 = [];
               end
               
               %Point is on an edge in x-space
               if abs(particle(index,2) - (RE - hx))<=10^(-12) || abs(particle(index,2)-LE) <= 10^(-12)
                   
                   index_wrapped2 = find(abs(particle(:,2) - particle(index,2))- Lx<= 10^(-12) ...
                       & abs(particle(:,2) - particle(index,2))-(Lx - hx)>= -10^(-12) ...
                       &  abs(particle(:,3) - particle(index,3))-hv<= 10^(-12)...
                       &  abs(particle(:,3) - particle(index,3))>= 0 );
               
               else 
                   
                   index_wrapped2 = [];
                   
               end
               
               %Point is on a corner
               if abs(particle(index,2) - (RE - hx))<=10^(-12) || abs(particle(index,2)-LE) <= 10^(-12) &&...
                       abs(particle(index,3) - (TE - hv))<= 10^(-12) || abs(particle(index,3) - BE)<=10^(-12)
                   
                   index_wrapped3 = find(abs(particle(:,2) - particle(index,2))- Lx<= 10^(-12) ...
                       & abs(particle(:,2) - particle(index,2))-(Lx - hx)>= -10^(-12) ...
                       &  abs(particle(:,3) - particle(index,3)) - Lv<= 10^(-12)...
                       &  abs(particle(:,3) - particle(index,3))-(Lv-hv)>=-10^(-12) );
               
               else 
                   
                   index_wrapped3 = [];
                   
               end
               
               index_neighbors = [index_neighbors; index_wrapped1; index_wrapped2; index_wrapped3];
               
           end
           
           index_neighbors(index_neighbors == index) = [];
           
           C_neighbors = particle(index_neighbors,1);
           
           C_neighbors_total = sum(C_neighbors(C_neighbors>0));
           
           if C_neighbors_total == 0 
               C_neighbors_total = Inf;
           end
           
           C_neighbors(C_neighbors<0 ) = 0;
           
           df = particle(negative_cell_index(j),1);
           
           Df = C_neighbors/C_neighbors_total*df; 
           
           particles_positive(index_neighbors) = particles_positive(index_neighbors)+Df;
           
           particles_positive(index, 1) = particles_positive(index, 1) -sum(Df);  
           
       end
       
       
       particle(:,1) = particles_positive;

       n = n+1;
       

       if n>=10
            particle(abs(particle(:,1))<10^(-10), 1) = 0;
       end
       
       if n >= 105
        error('too many iterations of positivity preservation')
       end
       
       
   end
       
       
       

end