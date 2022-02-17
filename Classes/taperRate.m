classdef taperRate
    
    methods(Static) 

        %==================================================================
        %>@brief    This construction is based on the width of dendritic segments at bifurcations, taken as a function of path
        %>          distance to the soma. Dendritic tapering is a measure of the change in width along a dendritic 
        %>          segment from node to node.
        %>@param none
        %>@retval generates a dendrogram and a distance matrix
        %==================================================================  
         function get_taper_rate()
            load(fullfile(pwd, CreateTree.save_directory('trees_vector.mat')));
            taper_rate_dist_matrix= taperRate.get_taper_rate_dist_matrix(trees_vector)
            save(fullfile( pwd , CreateTree.save_directory('taper_rate_dist_matrix')),'taper_rate_dist_matrix' );
            get_dendrogram(taper_rate_dist_matrix,'Taper_Rate');
         end 
         
        function taper_rate_dist_matrix= get_taper_rate_dist_matrix(trees_vector)
            vect_taper_rate =cell(1,length(trees_vector));
            vect_taper_rate{1,length(trees_vector)} = [];
            area_taper_rate = zeros(length(trees_vector),2);
            for i=1:length(trees_vector)          
               disp(['processing Taper Rate for tree :   ,' num2str(i), ' out of ' ,num2str(length(trees_vector))]);
                tree=normalize_tree(trees_vector(i));          
                 trv = taperRate.get_taper_rate_sholl_descriptor(tree);
                 vect_taper_rate{i}=trv;                
            end 
            save(fullfile( pwd , CreateTree.save_directory('vect_taper_rate')),'vect_taper_rate' );
            assignin('base','vect_taper_rate',vect_taper_rate);
            taper_rate_dist_matrix = zeros(length(trees_vector),length(trees_vector));
            for i=1:length(trees_vector)
                 for j=(i+1):length(trees_vector)
                   dtaperrate=DM.get_distance_of_plots_grid_version( vect_taper_rate{i}(:,1),vect_taper_rate{i}(:,2),vect_taper_rate{j}(:,1),vect_taper_rate{j}(:,2) ); 
                   taper_rate_dist_matrix(i,j) = dtaperrate;
                   taper_rate_dist_matrix(j,i) = dtaperrate;
                 end 
                 area_taper_rate(i,1) = trapz(vect_taper_rate{i}(:,1),vect_taper_rate{i}(:,2));
                if vect_taper_rate{i}(1,1) == 0 && vect_taper_rate{i}(1,2) ~= 0
                    area_taper_rate(i,2) = vect_taper_rate{i}(1,2);
                else
                    area_taper_rate(i,2) = vect_taper_rate{i}(end,2);               
                end 
                 [m , b ] = get_linear_regression(vect_taper_rate{i}(:,1),vect_taper_rate{i}(:,2));
                area_taper_rate(i,3) = m;
                area_taper_rate(i,4) = b;  
                
            end  
            save(fullfile( pwd , CreateTree.save_directory('area_taper_rate')),'area_taper_rate' );
        end

        function taper_rate_vector=get_taper_rate_sholl_descriptor(tree)
         
          nodesVector = tree.get_nodesVector();
          % get branches and leafs
          predicate_all = ones(1,length(tree.nodesVector));
          predicate1 = tree.is_branching_point();
          predicate2 = tree.is_leaf();
          predicate= predicate1 + predicate2;
          % get all the radii i.e the distances between the soma to all
          % bifricating points
            rad=[];
            root_coord = tree.get_root().coordinates;
            intrinsic_matrix = Spectral_descriptors.intrinsic_distance_matrix( tree , predicate_all);
                for i =1 :length(predicate)
                  if predicate(i)==1
                        taper_rate_vector(i,1) = intrinsic_matrix(1,i);
                        taper_rate_vector(i,2) = nodesVector(i).radius;
                  end 
                end  
%           delete zero rows
           taper_rate_vector( ~any(taper_rate_vector,2), : ) = [];
           taper_rate_vector= sortrows(taper_rate_vector,1); 
           maxval = max(taper_rate_vector(:,1));
           taper_rate_vector(:,1) = taper_rate_vector(:,1)/maxval;
        end
        
    end
end
