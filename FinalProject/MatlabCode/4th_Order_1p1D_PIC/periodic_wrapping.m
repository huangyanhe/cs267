%Need to rethink how I am dealing with periodicity Maybe I should write
%a function that wraps periodically around the boundary based on the
%number of points you need around the boundary and the dimension that
%you are in.

function particle_wrapped = periodic_wrapping(particle_one_domain, h_vector, domain_specs, size_of_extension)

 indices = domain_specs(:,4);

for j = 1:length(indices)
    
    Index = indices(j); 
    
    
    h = h_vector(j);

        length_of_domain = abs(domain_specs(j,2) - domain_specs(j,1));
    
        particle_extra_large = particle_one_domain(particle_one_domain(:, Index) >= domain_specs(j, 2)-h*size_of_extension ,:);
        particle_extra_small = particle_one_domain(particle_one_domain(:, Index) <= domain_specs(j, 1) + h*size_of_extension ,:);
        particle_extra_large(:, Index) = particle_extra_large(:, Index) - length_of_domain;
        particle_extra_small(:, Index) = particle_extra_small(:, Index) + length_of_domain;
        
        particle_one_domain = [particle_one_domain; particle_extra_small; particle_extra_large];
        
    
end


particle_wrapped = particle_one_domain;

end
