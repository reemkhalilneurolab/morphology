classdef Tmd_classical
          
    methods(Static)

        function get_tmd_dendrogram()
            load(fullfile(pwd, CreateTree.save_directory('trees_vector.mat')));
            tmd_dist_matrix_classical= Tmd_classical.get_tmd_dist_matrix(trees_vector);
            assignin('base','tmd_dist_matrix_classical',tmd_dist_matrix_classical);
            save(fullfile( pwd , CreateTree.save_directory('tmd_dist_matrix_classical')),'tmd_dist_matrix_classical' );
%             DM.get_clusters(tmd_dist_matrix);
            get_dendrogram(tmd_dist_matrix_classical,'tmd_classical');
        end 
        
        function tmd_dist_matrix= get_tmd_dist_matrix(trees_vector)
            TMD_diagram_classical = cell(1,length(trees_vector)); 
            persistence_pairs = cell(1,length(trees_vector)); 
            for i=1:length(trees_vector)
                mytree1=trees_vector(i);
               mytree=normalize_tree(mytree1);
               disp(['TMD_classical diagram for tree :   ,' num2str(i)]);

               diameter=mytree.get_tree_diameter();
               parts = 10;
               rad = zeros(parts,1);
               for j =1:parts
                  rad(j,1) = j*(diameter/parts);
               end  
               rad = sort(rad);
                TMD_diagram_classical{i} =  Tmd_classical.TMD_classical(mytree);
            end 
            save(fullfile( pwd , CreateTree.save_directory('TMD_diagram_classical')),'TMD_diagram_classical' );
            tmd_dist_matrix = zeros(length(trees_vector),length(trees_vector));
            for i=1:length(trees_vector)
                l1 = barcodeToLandscape(TMD_diagram_classical{i});
                 for j=(i+1):length(trees_vector)
                    disp(['calculating persistance distance between  :   ' num2str(i), ' and ' num2str(j)]);
                        l2 = barcodeToLandscape(TMD_diagram_classical{j});
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function  TMD_diagram = TMD_classical(mytree)
          
            %First we need to get a distance from the soma to all the nodes:
            distnace_function = zeros(1,mytree.get_num_nodes());
            root = mytree.get_root();

            disp(mytree.get_num_nodes())

            for i = 1:mytree.get_num_nodes()
                node = mytree.Get_Node(i);
                distnace_function(i) = sqrt( ( node.coordinates(1)-root.coordinates(1) )^2 + ( node.coordinates(2)-root.coordinates(2) )^2 + ( node.coordinates(3)-root.coordinates(3) )^2 );
            end   

            %Now we want to get branch decomposition using TMD_branch_decomposition
            pairings = Tmd_classical.TMD_branch_decomposition( mytree,distnace_function );
            %And now, for every pairing, we want to check get the values of
            %function:
                       
            TMD_diagram = [];
            
            for i=1:length(pairings)
                if ( pairings(i,2) ~= -1 )
                    TMD_diagram = [ TMD_diagram ; [distnace_function(pairings(i,1)) , distnace_function(pairings(i,2))] ];   
                else
                    TMD_diagram = [ TMD_diagram ; [distnace_function(pairings(i,1)) , 0] ];
                end
            end
            TMD_diagram;
        end
        
    end
end

