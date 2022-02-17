classdef Tmd_sholl
          
    methods(Static)        

        function get_tmd_dendrogram()
            load(fullfile(pwd, CreateTree.save_directory('trees_vector.mat')));
            tmd_dist_matrix_sholl= Tmd_sholl.get_tmd_dist_matrix(trees_vector);
            assignin('base','tmd_dist_matrix_sholl',tmd_dist_matrix_sholl);
            save(fullfile( pwd , CreateTree.save_directory('tmd_dist_matrix_sholl')),'tmd_dist_matrix_sholl' );
            get_dendrogram(tmd_dist_matrix_sholl,'tmd_sholl');
        end 
        
        function tmd_dist_matrix= get_tmd_dist_matrix(trees_vector)
            TMD_diagram_sholl = cell(1,length(trees_vector)); 
            persistence_pairs = cell(1,length(trees_vector)); 
            for i=1:length(trees_vector)
                mytree1=trees_vector(i);
               mytree=normalize_tree(mytree1);
               disp(['TMD_sholl diagram for tree :   ,' num2str(i)]);

               diameter=mytree.get_tree_diameter();
               parts = 10;
               rad = zeros(parts,1);
               for j =1:parts
                  rad(j,1) = j*(diameter/parts);
               end  
               rad = sort(rad);
               TMD_diagram_sholl{i} =  Tmd_sholl.sholl_tmd(rad, mytree);
            end 
            save(fullfile( pwd , CreateTree.save_directory('TMD_diagram_sholl')),'TMD_diagram_sholl' );
            tmd_dist_matrix = zeros(length(trees_vector),length(trees_vector));
            for i=1:length(trees_vector)
                l1 = barcodeToLandscape(TMD_diagram_sholl{i});
                 for j=(i+1):length(trees_vector)
                    disp(['calculating persistance distance between  :   ' num2str(i), ' and ' num2str(j)]);
                        l2 = barcodeToLandscape(TMD_diagram_sholl{j});
                        dist = landscapeDistance(l1, l2, 2);
                    tmd_dist_matrix(i,j) = dist;
                    tmd_dist_matrix(j,i) = dist;
                 end 
            end        
        end 
           
          
            function  branch_decomposition = TMD_branch_decomposition(mytree,function_on_vertices)
                %First we find all the leafs:
                is_leaf_ = mytree.is_leaf();

                %First, put all the leafs into the queue
                branches_function = [];
                branches_id = [];
                nodes = [];  
                for i = 1:length(is_leaf_)
                    if ( is_leaf_(i) == 1 )            
                        nodes = [ nodes mytree.Get_Node(i) ];
                        branches_function(i) = function_on_vertices(i);
                        branches_id(i) = i;
                    else
                        branches_function(i) = Inf;
                        branches_id(i) = -1;
                    end
                end

                branch_decomposition = [];


                %Now, starting from leafs, move up:
                while ( ~isempty(nodes) )      
                    new_nodes = [];

                    %         disp( "Here are the nodes we consider." )
                    %         for i = 1:length( nodes )
                    %             disp( nodes(i).id );
                    %         end

                    for i = 1:length( nodes )
                        node= nodes(i);                            

                        if ( node.parent == -1 )
                            %We have reached the root and we are done. Break;
                            branch_decomposition = [ branch_decomposition; branches_id( node.id ), -1 ];
                            break;
                        end

                        %find a parent of the node:
                        parent = mytree.Get_Node(node.parent);

                        if ( branches_function(parent.id) ~= Inf )
                            %The parent of this node has already been assigned, continue. 
                            %disp("Parent has already been asigned.")
                            continue;
                        end

                        %check if all the children have assigned the value of the function:
                        children = mytree.get_children(parent);
                        all_assigned = 1;                     

                        max_child = children(1);
                        for i=1:length( children )
                            if ( branches_function( children(i).id ) == Inf )
                                all_assigned = 0;
                            else
                                if ( branches_function( max_child.id ) < branches_function( children(i).id ) )
                                    max_child = children(i);
                                end
                            end
                        end

                        if ( all_assigned == 1 )
                            %all the children have asigned value of the funciton. Now, one
                            %branch will have to continue, and the rest will terminate:
                            for i=1:length( children )
                                if ( children(i) == max_child )
                                    %This branch will continuie
                                    branches_function( parent.id ) = branches_function( max_child.id );
                                    branches_id( parent.id ) = branches_id( max_child.id );
                                    %We add the node to the queue                        
                                    new_nodes = [new_nodes, parent];
                                else 
                                    %We kill this branch, i.e. it is being added to the
                                    %branch decomposition
                                    branch_decomposition = [ branch_decomposition; branches_id( children(i).id ), parent.id ];
                                end
                            end
                        else
                            %Do nothing, wait till all the children are assigned. 
                        end
                    end
                    nodes = unique(new_nodes);
                end
         end

%%%%%%%%%%%%%%%%%%%%%%
       
        function persistence_pairs = sholl_tmd(rad, tree)
          persistence_pairs =[];
          nodesVector = tree.get_nodesVector();
          pridicate_all = ones(1,length(tree.nodesVector));
          root = tree.get_root();
          for i = 1:tree.get_num_nodes()
                node = tree.Get_Node(i);
                distnace_function(i) = sqrt( ( node.coordinates(1)-root.coordinates(1) )^2 + ( node.coordinates(2)-root.coordinates(2) )^2 + ( node.coordinates(3)-root.coordinates(3) )^2 );
          end   
          pairings = Tmd_sholl.TMD_branch_decomposition( tree,distnace_function );
          root_coord = [root.coordinates(1) root.coordinates(2) root.coordinates(3)];          
          euclid_matrix = Spectral_descriptors.Euclidean_distance_matrix(tree,pridicate_all);
       for f=1:length(rad)
          for i =1:size(pairings,1) 
              if pairings(i,2) == -1             
                   dist = 0;
              else
                  dist = distance_of_vectors( nodesVector(pairings(i,2)).coordinates , root_coord ) ;
              end 
             if dist < rad(f) 
                ending = pairings(i,1);
                this_branch = [];
                while ( ending ~= pairings(i,2) )  
                    this_branch = [ this_branch , ending ];
                    ending = nodesVector(ending).parent;                    
                end 
                if ending ~=-1                    
                   this_branch = [ this_branch , ending ];
                end 
                cnt = 0;
                for j=1:length(this_branch)
                   if euclid_matrix(1,this_branch(1,j)) > rad(f)
                      if rad(f) < dist 
                      persistence_pairs = [persistence_pairs ;[rad(f) ,dist] ] ;
                      else
                      persistence_pairs = [persistence_pairs ;[ dist,rad(f)] ] ;
                      end 
                      cnt = 1;
                      break;
                   end 
                end
                if cnt ==0 
                    if ending ~=-1   
                        if euclid_matrix(1,pairings(i,1)) < dist
                        persistence_pairs = [persistence_pairs ;[ euclid_matrix(1,pairings(i,1)) , dist ] ] ;
                        else
                        persistence_pairs = [persistence_pairs ;[ dist , euclid_matrix(1,pairings(i,1)) ] ] ;
                        end 
                    else 
                        if euclid_matrix(1,pairings(i,1)) < -1 
                        persistence_pairs = [persistence_pairs ;[ euclid_matrix(1,pairings(i,1)) , -1 ] ] ;
                        else
                        persistence_pairs = [persistence_pairs ;[-1, euclid_matrix(1,pairings(i,1)) ] ] ;
                        end 
                    end 
                end 
%                 disp(['pairs(i,1)   ,' num2str(pairings(i,1) ), ' pairs(i,2) ' ,num2str(pairings(i,2)), ' this_branch ', num2str(this_branch) ]);
             end 
           end             
        end 
     end   
        
    end
end

