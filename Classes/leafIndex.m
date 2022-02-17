classdef leafIndex
    
    methods(Static) 
        %==================================================================
        %>@brief   This construction counts the number of terminations emanating from every node, plotted as
        %>         a function of radial distance. If the node is taken to be the soma, this number is the
        %>         total number of leaves, if the node is taken to be a leaf, the value is one (that leaf).
        %>@param none
        %>@retval generates a dendrogram and a distance matrix
        %==================================================================  
        
         function get_leaf_index()

            load(fullfile(pwd, CreateTree.save_directory('trees_vector.mat')));
            leaf_index_dist_matrix= leafIndex.get_leaf_index_dist_matrix(trees_vector);
            save(fullfile( pwd , CreateTree.save_directory('leaf_index_dist_matrix')),'leaf_index_dist_matrix' );
            get_dendrogram(leaf_index_dist_matrix,'leaf_index');
         end 
         
        function leaf_index_dist_matrix= get_leaf_index_dist_matrix(trees_vector)
            area_leaf_index = zeros(length(trees_vector),2);
            stats =[];
            vect_leaf_index =cell(1,length(trees_vector));
            vect_leaf_index{1,length(trees_vector)} = [];
            for i=1:length(trees_vector)
           
               disp(['processing leaf index for tree :   ,' num2str(i), ' out of ' ,num2str(length(trees_vector))]); 
                %normalize the tree
                 tree=normalize_tree(trees_vector(i));
 
                 liv = leafIndex.get_leaf_index_sholl_descriptor(tree);
                 med =  median( liv(:,6)) ; %median of raduis/leafindex
                 avrg = mean(liv(:,6)) ; % mean of raduis/leafindex
                 stddev= std(liv(:,6));
                 stats = [stats; med avrg stddev] ; 
                 vect_leaf_index{i}=liv(:,[4,2]);  
            end 
            assignin('base','vect_leaf_index',vect_leaf_index); 
            save(fullfile( pwd , CreateTree.save_directory('vect_leaf_index')),'vect_leaf_index' );
            leaf_index_dist_matrix = zeros(length(trees_vector),length(trees_vector));
            for i=1:length(trees_vector)
                 for j=(i+1):length(trees_vector)
                   dLeafIndex=DM.get_distance_of_plots_grid_version( vect_leaf_index{i}(:,1),vect_leaf_index{i}(:,2),vect_leaf_index{j}(:,1),vect_leaf_index{j}(:,2) ); 
                   leaf_index_dist_matrix(i,j) = dLeafIndex;
                   leaf_index_dist_matrix(j,i) = dLeafIndex;
                 end 
                area_leaf_index(i,1) = trapz(vect_leaf_index{i}(:,1),vect_leaf_index{i}(:,2));
                if vect_leaf_index{i}(1,1) == 0 && vect_leaf_index{i}(1,2) ~= 0
                    area_leaf_index(i,2) = vect_leaf_index{i}(1,2);
                else
                    area_leaf_index(i,2) = vect_leaf_index{i}(end,2);               
                end  
                
               [m , b ] = get_linear_regression(vect_leaf_index{i}(:,1),vect_leaf_index{i}(:,2));
                area_leaf_index(i,3) = m;
                area_leaf_index(i,4) = b;  
                 
            end   
            save(fullfile( pwd , CreateTree.save_directory('area_leaf_index')),'area_leaf_index' );

        end
        
        function all_path_vector=get_leaf_index_sholl_descriptor(tree)
         
            all_path= tree.get_paths_from_leafs_to_root();
            len = length(all_path);
            mypath=[];
            for i=1:len
              mypath = [mypath all_path{i}] ;   
            end 
            assignin('base','mypath',mypath); 
            [num_of_leafs,nodes_id]=hist(mypath,unique(mypath));
            
            % node id ---num of leafs ---bif and leaf flag ---distance from
            % node to soma---node raduis ---node raduis/leafindex
            all_path_vector = [nodes_id; num_of_leafs]; 
            all_path_vector =all_path_vector';
            [r,c]=size(all_path_vector);
            root_coord = tree.get_root().coordinates;
            for i=1:r
                    node = tree.Get_Node(all_path_vector(i,1));
                    children = tree.get_children(node);
                    if ( length(children) > 1 || length(children) == 0 )
                        all_path_vector(i,3) = 1;
                        all_path_vector(i,4) = distance_of_vectors( node.coordinates , root_coord );
                        all_path_vector(i,5) = node.radius; 
                        all_path_vector(i,6) = all_path_vector(i,5)/ all_path_vector(i,2);
                    end
            end
%          all_path_vector : c1 = node id ----c2 = leaf index ---c3 = isleaf, is branch ---c4=raduis
            delete_rows = all_path_vector(:,3) == 0;   % for which rows is column 3 equal to zero
            all_path_vector(delete_rows,:) = [];   % remove those rows
            all_path_vector= sortrows(all_path_vector,4); %  sort by raduis  

            
                
        end
        
    end
end

