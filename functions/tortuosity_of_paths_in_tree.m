%==========================================================================
%> @brief  tortuosity_of_paths_in_tree return a vector of values if
%tortuosity for every path in the paths cell array. To get the paths,
%please use get_paths_from_leafs_to_root method of the Tree class. 
%>
%> @retval vector of tortuosity of paths in the paths input variable
% =========================================================================
function tortuosity_vector = tortuosity_of_paths_in_tree( tree , paths )
    tortuosity_vector = zeros( 1,length( paths ) );
    length_of_path_vector = zeros( 1,length( paths ) );
    %for every path:
    for i = 1:length( paths )
        %compute a distance between the endpoints of path:
        path = paths{i}
        
        %First, let us get the first and the last node, and compute the
        %distnace between them:
        first_node = tree.Get_Node( path(1) );
        last_node = tree.Get_Node( path( length(path)  ) );
        distance_between_endpoints = distance_of_vectors( first_node.coordinates , last_node.coordinates )
        
        %Now compute the length of a path, that is a sum of lengths of all
        %the intermediate line segemnts:
        length_of_path = 0;
        for j = 2:length( path )
            path(j-1)
            path(j)
            previous = tree.Get_Node( path(j-1) );
            current =  tree.Get_Node( path(j) );
            length_of_path = length_of_path + distance_of_vectors( previous.coordinates , current.coordinates );
        end
        %Now we have all to compute tortuosity:
%         tortuosity_function=[length_of_path/distance_between_endpoints ; length_of_path];
%         tortuosity_vector = [ tortuosity_vector,tortuosity_function  ]; 
        tortuosity_vector(i) = length_of_path/distance_between_endpoints ;
        length_of_path_vector(i)=length_of_path;
        
    end
    myvector=[tortuosity_vector;length_of_path_vector]
    myvector=myvector';     
    tortuosity_vector=sortrows(myvector);
   
end

