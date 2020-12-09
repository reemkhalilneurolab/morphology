classdef Spectral_descriptors
    
    methods(Static) 
      
%         ==================================================================
        %>@brief   
        %>@param 
        %>@retval
        %==================================================================                     
        function matrix = Euclidean_distance_matrix(tree,predicate)
             matrix = zeros( length(tree.nodesVector) , length(tree.nodesVector));
%             matrix =sparse( length(tree.nodesVector) , length(tree.nodesVector));
                for i = 1:length(tree.nodesVector)     
                    if ( predicate(i) == 1 )
                        for j = i+1:length(tree.nodesVector)
                            if ( predicate(i) == 1 )
                                d = distance_of_vectors( tree.nodesVector(i).coordinates , tree.nodesVector(j).coordinates );
                                matrix( i,j ) = d;
                                matrix( j,i ) = d;
                            end
                        end
                    end
                end 
            matrix( ~any(matrix,2), : ) = [];  %remove zero rows
            matrix( :, ~any(matrix,1) ) = [];  %remove zero columns
            colums_to_consider = zeros( 1 , sum(predicate) );
            nr = 1;
                for i=1:length(predicate)
                    if ( predicate(i) == 1 )
                        colums_to_consider(nr) = i;
                        nr = nr + 1;
                    end
                end
            matrix = matrix( colums_to_consider,colums_to_consider );
        end

        %==================================================================
        %>@brief   
        %>@param 
        %>@retval
        %==================================================================
        function eigenvalues = Euclidean_distance_between_nodes(tree,predicate)
            matrix = Euclidean_distance_matrix(tree,predicate);
            eigenvalues = eig(matrix);
        end

        %==================================================================
        %>@brief this should be invariant with respect to the size of the tree. 
        %>@param 
        %>@retval
        %==================================================================
        function eigenvalues_normalized = normalized_spectral_radius_decay_profile_Euclidean_distance(tree)
            predicate = ones(1,length(tree.nodesVector));
            eigenvalues = Euclidean_distance_between_nodes(tree,predicate);
            eigenvalues = sort( abs(eigenvalues) );
            eigenvalues_normalized = normalize_vector_01( eigenvalues );
        end

        %==================================================================
        %>@brief   
        %>@param 
        %>@retval
        %==================================================================        
        function spec_rad = spectral_radius_of_tree_Euclidean_distance(tree)
            predicate = ones(1,length(tree.nodesVector));
            eigen = Spectral_descriptors.Euclidean_distance_between_nodes(tree,predicate);
            spec_rad = max(abs(eigen));
        end

        %==================================================================
        %>@brief   
        %>@param 
        %>@retval
        %==================================================================
        function spec_rad = spectral_radius_of_tree_Euclidean_distance_branches_only(tree)
            predicate = tree.is_branching_point();
            eigen = Spectral_descriptors.Euclidean_distance_between_nodes(tree,predicate);
            spec_rad = max(abs(eigen));
        end  

        %==================================================================
        %>@brief   
        %>@param 
        %>@retval
        %==================================================================
        function spec_rad = spectral_radius_of_tree_Euclidean_distance_leafs_only(tree)
            predicate = tree.is_leaf();
            eigen = Spectral_descriptors.Euclidean_distance_between_nodes(tree,predicate);
            spec_rad = max(abs(eigen));
        end
        
     
        function d = intrinsic_distance_matrix( mytree , predicate )
            G = GRAPH.convert_to_matlab_graph(mytree);            
            d = distances(G);
            %Now, we will remove the columns corresponding to the
            %predicates. I assume here that enumeration of nodes in the
            %graph and the enumeration of columns in the matrix are the
            %same. This should be checked!!!
            %Checked ..all good 
            colums_to_consider = zeros( 1 , sum(predicate) );
            nr = 1;
                for i=1:length(predicate)
                    if ( predicate(i) == 1 )
                        colums_to_consider(nr) = i;
                        nr = nr + 1;
                    end
                end
            d = d( colums_to_consider,colums_to_consider );
        end

        function eigenvalues = intrinsic_distance_between_nodes(mytree,predicate)
            d = Spectral_descriptors.intrinsic_distance_matrix( mytree , predicate );
            eigenvalues = eig(d);
        end
        
      
        function spec_rad = spectral_radius_of_tree_intrinsic_distance(tree)
            predicate = ones(1,length(tree.nodesVector));
            eigen = Spectral_descriptors.intrinsic_distance_between_nodes(tree,predicate);
            spec_rad = max(abs(eigen));
        end
        

        function spec_rad = spectral_radius_of_tree_intrinsic_distance_branches_only(tree)
            predicate = tree.is_branching_point();
            eigen = Spectral_descriptors.intrinsic_distance_between_nodes(tree,predicate);
            spec_rad = max(abs(eigen));
        end
        
        function spec_rad = spectral_radius_of_tree_intrinsic_distance_leafs_only(tree)
            predicate = tree.is_leaf();
            eigen = Spectral_descriptors.intrinsic_distance_between_nodes(tree,predicate);
            spec_rad = max(abs(eigen));
        end

        function spec_rad = spectral_radius_of_tree_tortuosity_branches_only(tree)
            predicate = tree.is_branching_point();
            eigen = Spectral_descriptors.eigenvalues_of_toruosity_matrix(tree,predicate);
            spec_rad = max(abs(eigen));
        end

        function spec_rad = spectral_radius_of_tree_tortuosity_leafs_only(tree)
            predicate = tree.is_leaf();
            eigen = Spectral_descriptors.eigenvalues_of_toruosity_matrix(tree,predicate);
            spec_rad = max(abs(eigen));
        end            
    end
end

