classdef branchingPattern
   
    methods(Static)  
        %==================================================================
        %>@brief   This is an integer valued descriptor given by the number of bifurcations from which
        %>         we  subtract the number of leaves, all within a given radius r from the soma. 
        %>@param none
        %>@retval generates a dendrogram and a distance matrix
        %==================================================================        
        function get_branching_pattern()
            %load the trees vector
            load(fullfile(pwd, CreateTree.save_directory('trees_vector.mat')));
            %get the distance matrix
            branching_dist_matrix=branchingPattern.get_dist_matrix(trees_vector);    
            save(fullfile( pwd , CreateTree.save_directory('branching_dist_matrix')),'branching_dist_matrix' );
            %generate dendrogram
            get_dendrogram(branching_dist_matrix,'Branching');
        end 
         
        %==================================================================
        %>@brief  get_dist_matrix takes a cell array trees , normalizes
        %each of them and converts to a function. Distances between
        %fucntions is then computed and a distance matrix is returned. 
        %>@param a vector of trees
        %>@retval distance matrix of distances between fucntions
        %==================================================================

        function dist_matrix=get_dist_matrix(trees_vector)
            area_branching_pattern = zeros(length(trees_vector),2);
            vect_branching_pattern =cell(1,length(trees_vector));
            for i=1:length(trees_vector)
                disp(['processing branching distance for tree :   ,' num2str(i), ' out of ' ,num2str(length(trees_vector)), ' --started at ', datestr(now,'HH:MM:SS.FFF')]);
                %normalize the tree
                tree=normalize_tree(trees_vector(i));
                %tree.branches_and_leafs_distances_plot takes any function on the nodes and plot the values
                %of the following function f. f : r -> number of branching points
                %in the ball centered in the root and radius r minus number of
                [x1,y1] = tree.branches_and_leafs_distances_plot( tree.get_distances() );                 
                vect_branching_pattern{i}=[x1',y1'];            
            end 

            save(fullfile( pwd , CreateTree.save_directory('vect_branching_pattern')),'vect_branching_pattern' );
            dist_matrix = zeros(length(trees_vector),length(trees_vector));
            for i=1:length(trees_vector)
                 for j=(i+1):length(trees_vector)
                   %get distances between functions
                   d=DM.get_distance_of_plots_grid_version( vect_branching_pattern{i}(:,1),vect_branching_pattern{i}(:,2),vect_branching_pattern{j}(:,1),vect_branching_pattern{j}(:,2) );
                   dist_matrix(i,j) = d;
                   dist_matrix(j,i) = d;
                 end 
                %get area under the curve for vectorization
                area_branching_pattern(i,1) = trapz(vect_branching_pattern{i}(:,1),vect_branching_pattern{i}(:,2));
%                 phi(0) or phi(1)
                if vect_branching_pattern{i}(1,1) == 0 && vect_branching_pattern{i}(1,2) ~=0
                    area_branching_pattern(i,2) = vect_branching_pattern{i}(1,2);
                else
                    area_branching_pattern(i,2) = vect_branching_pattern{i}(end,2);               
                end
                
                [m , b ] = get_linear_regression(vect_branching_pattern{i}(:,1),vect_branching_pattern{i}(:,2));
                area_branching_pattern(i,3) = m;
                area_branching_pattern(i,4) = b;     
                
            end         
        save(fullfile( pwd , CreateTree.save_directory('area_branching_pattern')),'area_branching_pattern' );

        end
       
    end
end

