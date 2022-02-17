classdef tortuosity
    
    methods(Static)
        %==================================================================
        %>@brief  Tortuosity between two nodes is measured as the quotient 
        %>        of the path distance by the Euclidean distance. The tortuosity 
        %>        descriptor measures the mean tortuosity between all adjacent 
        %>        nodes of part of the neuron that is within radius r from the soma.
        %>@param none
        %>@retval generates a dendrogram and a distance matrix
        %==================================================================  
        function get_tortuosity()
            load(fullfile(pwd, CreateTree.save_directory('trees_vector.mat')));
            tortuosity_dist_matrix = tortuosity.get_tort_dist_matrix(trees_vector);
            get_dendrogram(tortuosity_dist_matrix,'Tortuosity');
        end 
        
        function tortuosity_dist_matrix = get_tort_dist_matrix(trees_vector)
            vect_tortuosity = cell(1,length(trees_vector));
            vect_tortuosity{1,length(trees_vector)} = [];            
            area_tortuosity = zeros(length(trees_vector),2);

            for i=1:length(trees_vector)
               disp(['processing Tortuosity for tree :   ,' num2str(i), ' out of ' ,num2str(length(trees_vector)), ' --started at ', datestr(now,'HH:MM:SS.FFF')]);
                 tree=normalize_tree(trees_vector(i));
                 tv = tortuosity.get_tort_func_r(tree);
                 vect_tortuosity{i}=tv;                
            end 
            
            save(fullfile( pwd , CreateTree.save_directory('vect_tortuosity')),'vect_tortuosity' );           
            assignin('base','vect_tortuosity',vect_tortuosity);            
            
            tortuosity_dist_matrix = zeros(length(trees_vector),length(trees_vector));
            for i=1:length(trees_vector)
                 for j=(i+1):length(trees_vector)
                   dtortuosity=DM.get_distance_of_plots_grid_version( vect_tortuosity{i}(:,1),vect_tortuosity{i}(:,2),vect_tortuosity{j}(:,1),vect_tortuosity{j}(:,2) ); 
                   tortuosity_dist_matrix(i,j) = dtortuosity;
                   tortuosity_dist_matrix(j,i) = dtortuosity;
                 end
                 
                 area_tortuosity(i,1) = trapz(vect_tortuosity{i}(:,1),vect_tortuosity{i}(:,2));                 
                 
                if vect_tortuosity{i}(1,1) == 0 && vect_tortuosity{i}(1,2) ~= 0
                    area_tortuosity(i,2) = vect_tortuosity{i}(1,2);
                else
                    area_tortuosity(i,2) = vect_tortuosity{i}(end,2);               
                end   
                
                [m , b ] = get_linear_regression(vect_tortuosity{i}(:,1),vect_tortuosity{i}(:,2));
                area_tortuosity(i,3) = m;
                area_tortuosity(i,4) = b;                
           
            end                     
             
             save(fullfile( pwd , CreateTree.save_directory('area_tortuosity')),'area_tortuosity' );
             save(fullfile( pwd , CreateTree.save_directory('tortuosity_dist_matrix')),'tortuosity_dist_matrix' );
             
             get_dendrogram(tortuosity_dist_matrix,'Tortuosity');

            
        end
      

        function grandmtx = get_tort_func_r(tree)
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
            average_tortuosity =[];
            grandmtx =[];
            for i=1:length(rad) % loop through all the radii 
            newpath =[];
            pathmtx=[];
            ratio =[] ;
            euclid_dist = [];
            intrinsic_dist = [];
                rw=1;
                for r=1:sz %loop though all the paths
                    mat=all_path{r};
                    cnt =1;
                    newpath =[];
                    
                    % loop througH the path and get a newpath vector for
                    % the points insdie the ball or on the sphere.
                    for d=size(mat,2):-1:1 
                         if euclid_matrix(1,mat(1,d)) <= rad(i) % if point on the path is inside the ball or on the sphere   
                            newpath(1,cnt) = mat(1,d);
                            newpath(2,cnt)= predicate_branch(1,mat(1,d)) ;
                            cnt=cnt+1;
                         end
                    end
                    
                    if euclid_matrix(1,newpath(1,end)) <= rad(i)
                        newpath(2,end)=1;% last point is always a new termination point if it is inside the ball
                        newpath(2,1)=1;
                    end                    
                    
%                    get all the branching points on that path
                     ind = newpath(2,:) == 1;
                     A1 = newpath(1,ind);
                                          
                     if ~isempty(A1)
                       for h = 1:length(A1)-1
                           pathmtx(rw,:)= [A1(1,h) A1(1,h+1)];
                           rw=rw+1;
                       end 
                     end                  
                                         
                end 
                
                if ~isempty(pathmtx)
                    B = unique(pathmtx,'rows');
                    for e=1:size(B,1)
                        B(e,3) = intrinsic_matrix(B(e,1),B(e,2));
                        B(e,4) = euclid_matrix(B(e,1),B(e,2));
                        
                        if B(e,4)== 0
                            B(e,5) = 0;
                        else
                           B(e,5) = B(e,3)/B(e,4);  
                        end
                    end 
                    grandmtx = [grandmtx ; rad(i) sum(B(:,5))/size(B,1)];
  
                end                
                
            end       

         end     
    

  end
end

