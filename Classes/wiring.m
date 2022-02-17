classdef wiring
    
    methods(Static)
        %==================================================================
        %>@brief  This construction measures the total wiring of a neuron
        %>        within a given sphere of radius r centered at the origin. 
        %>        This is the sum of the path distances of all dendritic 
        %>        branches within that sphere.                
        %>@param none
        %>@retval generates a dendrogram and a distance matrix
        %================================================================== 
        function get_wiring()
            load(fullfile(pwd, CreateTree.save_directory('trees_vector.mat')));
            wiring_dist_matrix= wiring.get_wiring_dist_matrix(trees_vector);
            save(fullfile( pwd , CreateTree.save_directory('wiring_dist_matrix')),'wiring_dist_matrix' );
            get_dendrogram(wiring_dist_matrix,'wiring');
        end 


        function wiring_dist_matrix= get_wiring_dist_matrix(trees_vector)
            area_wiring = zeros(length(trees_vector),2);
            for i=1:length(trees_vector)
               disp(['processing wiring for tree :   ,' num2str(i), ' out of ' ,num2str(length(trees_vector)), ' --started at ', datestr(now,'HH:MM:SS.FFF')]);
               tree=normalize_tree(trees_vector(i));
               wv = wiring.get_wiring_func_r(tree);
               vect_wiring{i}=wv;
                
            end 

            save(fullfile( pwd , CreateTree.save_directory('vect_wiring')),'vect_wiring' );
            assignin('base','vect_wiring',vect_wiring);
            wiring_dist_matrix = zeros(length(trees_vector),length(trees_vector));
            for i=1:length(trees_vector)
                 for j=(i+1):length(trees_vector)
                   dwiring=DM.get_distance_of_plots_grid_version( vect_wiring{i}(:,1),vect_wiring{i}(:,2),vect_wiring{j}(:,1),vect_wiring{j}(:,2) ); 
                   wiring_dist_matrix(i,j) = dwiring;
                   wiring_dist_matrix(j,i) = dwiring;
                 end 
                 area_wiring(i,1) = trapz(vect_wiring{i}(:,1),vect_wiring{i}(:,2));
                if vect_wiring{i}(1,1) == 0 && vect_wiring{i}(1,2) ~= 0
                    area_wiring(i,2) = vect_wiring{i}(1,2);
                else
                    area_wiring(i,2) = vect_wiring{i}(end,2);               
                end
                
                [m , b ] = get_linear_regression(vect_wiring{i}(:,1),vect_wiring{i}(:,2));
                area_wiring(i,3) = m;
                area_wiring(i,4) = b;    
                 
                
            end 
             save(fullfile( pwd , CreateTree.save_directory('area_wiring')),'area_wiring' );
        end
      
        function wiring_vect = get_wiring_func_r(tree)

            all_path=tree.get_paths_from_leafs_to_root();
            nodesVector = tree.get_nodesVector();
            pridicate_all = ones(1,length(tree.nodesVector));
            predicate_branch = tree.is_branching_point();
            predicate_leaf = tree.is_leaf();
            predicate=predicate_branch + predicate_leaf;

            % get all the radii i.e the distances between the soma to all
            % bifricating points
            rad=[];
            root_coord = tree.get_root().coordinates;
                for i =1 :length(predicate)
                  if predicate(i)==1
                     rad = [rad; distance_of_vectors( nodesVector(i).coordinates , root_coord )];  
                  end 
                end  
                rad= sort(rad);        
            %get euclidean distance matrix between all points/nodes/markers
            euclid_matrix = Spectral_descriptors.Euclidean_distance_matrix(tree,pridicate_all);
            % get intrinsic distance matrix between all points/nodes/markers
            intrinsic_matrix = Spectral_descriptors.intrinsic_distance_matrix( tree , pridicate_all);
            mat=[];
            pair_euclid_dist=0;
            pair_intrinsic_dist=0;

            sz = length(all_path);
          
            wiring_vect =[];
            for i=1:length(rad) % loop through all the radii 
                ratio =[] ;
                euclid_dist = [];
                intrinsic_dist = [];
                total_wiring =0;
                pair_mat =[];           
            
                for r=1:sz %loop though all the paths
                    mat=all_path{r};
                    for j=length(mat):-1:1 % start at the soma ..the all_path vector starts at the leaf and the last index is the soma ..we start at the last index
                    dist =euclid_matrix(1,mat(1,j));
                        if dist > rad(i) || j == 1
                            break;
                        end 
                        if all_path{r}(1,j) ~= 1
                         pair_mat = [pair_mat ; [mat(1,j+1) mat(1,j) intrinsic_matrix(mat(1,j+1),mat(1,j)) ]];
                        end 
%                     
                    end                 
                end 
                if ~isempty(pair_mat)
                    out = unique(pair_mat,'rows'); % delete duplicate paths
                    total_wiring = sum(out(:,3));% add all the distances for the unique paths
                end 
                wiring_vect = [wiring_vect ; rad(i) total_wiring]; 
            end
        
         end     
    
  end
end

