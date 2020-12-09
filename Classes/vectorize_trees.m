classdef vectorize_trees    
    methods(Static)
         
        function  tortuosity=tortuosity(mytree)
            %get the flux vectors
            paths= mytree.get_paths_from_leafs_to_root();
            tor=TORTUOSITY.tortuosity_of_paths_in_tree( mytree , paths );
            tor=tor(:,1)';
            tortuosity = normalize_vector_01(tor);
        end
        
        
        %Add comments
        function  [n, s, l]=flux(mytree)
            %get the flux vectors
            [~ ,n, s, l]=FLUX.get_flux(mytree); 
            n=n(:,2)';
            s=s(:,2)';
            l=l(:,2)';
            n = normalize_vector_01(n);
            s = normalize_vector_01(s);
            l = normalize_vector_01(l);
        end
        
        %Add comment
        function  energy=energy( mytree, min_x , max_x , min_y , max_y , min_z , max_z , dx )
            energy = ENERGY.compute_energy_on_grid( mytree , min_x , max_x , min_y , max_y , min_z , max_z , dx );
            v = size(energy);
            energy = reshape( energy , [1,v(1)*v(2)*v(3)] );
            %As we want to have something rotational invariant, let us sort it.
            energy = sort(energy);
            energy = normalize_vector_01(energy);
        end
        
        %Add comment
        function  eigenvals=eigenvalues( mytree )                
            eigenvals = (GRAPH.eigenvalues_of_distance_matrix(mytree))';
            eigenvals = normalize_vector_01(eigenvals);
        end
        
        %Add comment
        function  distance=distances( mytree )        
           [distance,~] = mytree.branches_and_leafs_distances_plot( mytree.get_distances() );
           distance = normalize_vector_01( distance );
        end

        %Add comment
        function rescaled_vect = rescale_vector( vect , N )
           rescaled_vect=zeros(1,N);
           x=1;
           dx=length(vect)/N;  
           for i =1:N
               rescaled_vect(i)=vect(floor(x));
               x=x+dx;
           end
        end    
    end
end

