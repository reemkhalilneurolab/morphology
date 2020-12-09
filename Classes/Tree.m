
classdef Tree < matlab.mixin.Copyable
   
    properties
        root; %tree root
        numNodes; %total number of nodes
        searchNode; % parent node being searched by id so we can attach a child to it
        nodesVector; %a vector of all the nodes  
        pnodes;
        
        
    end
    
    methods (Access = public)
%%
        function this = Tree()
            %-------------------------------
            %constructor: sets root and search node to NaN and number of
            %nodes to zero
            %-----------------------------          
          this.root = NaN;
          this.searchNode=NaN;
          this.numNodes=0;
          this.nodesVector=[];
          pnodes=[];

        end
 %%            
        function Traverse(this, node)
            if isnan(node) 
            else 
                        this.pnodes=[this.pnodes node];
                        disp(node.id);
                        Traverse(this,node.leftChild);
                        Traverse(this,node.rightSibling);
            end
        end   
      
%%                      
        function traverseTree(this,rt)
                if isnan(rt)
                else
                    while ~isnan(rt)
                       disp(rt.id);
                            if ~isnan(rt.leftChild)
                                traverseTree(this,rt.leftChild);
                            end
                        rt = rt.rightSibling;
                    end 
               end 
        end
 %%   
 
 %%
         function q=add_nodes_to_queue (this, node,q)
                if isnan(node)                   
                else                    
                            q.add(node);
                            q=add_nodes_to_queue(this,node.leftChild,q);
                            q=add_nodes_to_queue(this,node.rightSibling,q);
                end                  
            end
 %%       
 %%
        function num_nodes = get_num_nodes(this)
            num_nodes = this.numNodes;
        end
 %%
 %%       
         function num_bif_nodes = get_num_bif_nodes(this)
             branches=this.is_branching_point();
            num_bif_nodes = sum(branches);
         end
 %%
 %%
        function num_leaf_nodes = get_num_leaf_nodes(this)
            leafs=this.is_leaf();
            num_leaf_nodes=sum(leafs);
        end
 %%
 %%
          function searchNode = get_searchNode(this)
            searchNode = this.searchNode;
          end
 %%
 %%
          function root = get_root(this)
            root = this.root;
          end
 %%
 %%         
           function nodesVector = get_nodesVector(this)
            nodesVector = this.nodesVector;
           end
 %%
 %% 
         function rmc=getRightmostChild(this,node)
             %----------------------------------
             %loops through all the children/siblilings of the left child             
             %----------------------------------
             current= node;
              while ~isnan(current.rightSibling)              
                current=current.rightSibling;
              end 
          rmc=current;
         end 
%%
%%   
        function addNode(this,node)    
               %-------------------------------
            %addnode: creates a node , sets node properties from vect
            %values. vect is a vector with 7 colums to keep it inline with
            %the swc file
            %vect=[Sample# StructureIdentifier	X Y	Z Radius parentsample]
            %addnode adds leftchild first attaches the rest of the children
            %as right siblings of the leftchild.
            %-----------------------------            
            %firstnode =root
            if this.numNodes==0
                 this.nodesVector = [this.nodesVector node];   
                 this.root = node;
                 this.root.parent=-1;
                 this.numNodes=this.numNodes+1;
                else
                  %every other node that is not a root  
                 this.searchNode=NaN;
                 %search for the node
                 searchParentNode(this,node.parent,this.root); 
                 %if found
                  if ( ~isnan(this.searchNode))
                          %check if the father has any children
                          %if the father has no children add the node as a left
                          %child
                         if isnan(this.searchNode.leftChild) 
                          this.searchNode.leftChild =node;
                          node.rightSibling = NaN;                  
                          %if the father has children check all the children and
                          %add the node as the right most node
                         else
                             rmc=getRightmostChild(this,this.searchNode.leftChild);                        
                             rmc.rightSibling =node;

                         end                     
                        this.nodesVector = [this.nodesVector node]; 
                        this.numNodes=this.numNodes+1;
                   else                 
                   end
            end                   
        end
%%
%%      
         function searchParentNode(this,newNodeParent,currentNode)
              %--------------------------------------------------
              %Recursively searches for the parent node and returns it so
              %we can add the node as one of the children
              %------------------------------------------------
              this.searchNode = this.nodesVector(newNodeParent);            
         end 
%%
%%        
         function nVector=getNodesVector(this)
               %----------------------------------------
               %returns a matrix of all the nodes               
               %----------------------------------
               nVector=this.nodesVector;
         end
%%
%%            
          function mynode = Get_Node(this, node_id) 
              if ( (node_id > 0) && (node_id <= this.numNodes) )
                  mynode = this.nodesVector(node_id) ;                
              else
                  mynode = nan;
              end
          end
%%
%%            
         function children=get_children(this,node)                
                children = [];                 
                temp_node = node.leftChild;                             
                while ~isnan(temp_node)
                   children = [children , temp_node];                        
                   temp_node=temp_node.rightSibling;
               end                
         end
%%
%%               
        function primary_trees = decompose_primary_branches(this)
                rt=this.get_root();
                children =this.get_children(rt); 
                primary_trees=[];
                for i=1:length(children)
                    if ~isnan(children(i).leftChild)
                        primary_trees = [primary_trees , tree];
                    end                                   
                end 
        end 
%%
%%        
       function is_branch = is_branching_point(this)
                is_branch = zeros(1,this.numNodes);
                for i=1:this.numNodes
                    node = this.Get_Node(i);
                    children = this.get_children(node);
                    if ( length(children)==2 && node.parent ~=-1)
                        is_branch(i) = 1;
                    end
                end
       end
%%
%%  
       function is_leaf = is_leaf( this )
                is_leaf = zeros(1,this.numNodes);
                for i=1:this.numNodes
                    node = this.Get_Node(i);
                    children = this.get_children(node);
                    if ( length(children) == 0 )
                        is_leaf(i) =1;
                    end
                end
       end
%%
%%
         %This function restricts the values obtained from get_distances 
          %method to all branching points.        
          function dist_bif_points = get_distance_bifurcation_points(this)
                max_dist=this.get_distances();
                is_branch = this.is_branching_point( );
                dist_bif_points = [];
                for i = 1:length( is_branch )
                    if ( is_branch(i) == 1 )
                        dist_bif_points = [dist_bif_points max_dist(i) ];
                    end
                end
          end
%%
%%
          %This function assign, for every node N of the tree, a value
          %equal to the maximal distance from the root to any point in the 
          %branch between the node N and the root. 
          function max_dist=get_distances(this)              
               %Load tree nodes in a queue
                nodes_q = this.add_nodes_to_queue(this.root,Queue('Node') ); 
                %create a q to process all the nodes and get distances
                dist_q = Queue('Node'); 
                %craete a zeros array to save all the distances
                max_dist = zeros(1,nodes_q.size());
                %add the root
                dist_q.add(this.root);
                while (~dist_q.isempty())
                    node= dist_q.pop();
                    children = this.get_children(node);
                    for i =1:length(children)
                         cdistance=sqrt((children(i).coordinates(1)-this.root.coordinates(1))^2+(children(i).coordinates(2)- this.root.coordinates(2))^2+(children(i).coordinates(3)-this.root.coordinates(3))^2);  
                             if cdistance > max_dist(node.id)
                                 max_dist(children(i).id)=cdistance;
                             else
                                 max_dist(children(i).id)=max_dist(node.id);
                             end
                         dist_q.add(children(i));  
                    end                  
                end              
          end
%%
%%
        %This function takes any given function on the nodes and return
        %its values on branching points and leafs. 
        function [distance_of_branches , distance_of_leafs] = branches_and_leafs_distances( this, function_on_nodes )
            is_branch = this.is_branching_point( );
            is_leaf = this.is_leaf( );            
            %now, for every branch and leaf we will get its distance from somma:
            distance_of_branches = []; %Please note that in this case I am treating root as a branch!!
            distance_of_leafs = [];
            root = this.root;
            for i=1:this.numNodes
                node = this.Get_Node(i);        
                if ( is_branch(i) == 1 )
                    distance_of_branches = [ distance_of_branches , function_on_nodes(i) ];
                end        
                if ( is_leaf(i) == 1 )
                    distance_of_leafs = [ distance_of_leafs , function_on_nodes(i) ];
                end        
            end
            distance_of_branches = sort( distance_of_branches );
            distance_of_leafs = sort( distance_of_leafs );
        end
%%
%%
        %This function takes any function on the nodes and plot the values
        %of the following function f. f : r -> number of branching points
        %in the ball centered in the root and radius r minus number of
        %leafs therein. 
        function [x,y] = branches_and_leafs_distances_plot( this , function_on_nodes )
            [distance_of_branches , distance_of_leafs] = this.branches_and_leafs_distances( function_on_nodes );
            x = [];
            y = [];
            value = 0;
            branches_iterator = 1;
            leafs_iterator = 1;
            while ( (branches_iterator <= length(distance_of_branches)) && (leafs_iterator <= length(distance_of_leafs)) )
                if ( distance_of_branches(branches_iterator) <= distance_of_leafs(leafs_iterator) )
                    %In this case a branch is encountered in a given radius. We
                    %increase the value.
                    x = [ x distance_of_branches(branches_iterator) ];
                    value = value + 1;          
                    branches_iterator = branches_iterator + 1;
                else
                    %In this case a leaf is encountered in a given radius. We
                    %decrease the value.
                    x = [ x , distance_of_leafs(leafs_iterator) ];
                    value = value - 1  ;
                    leafs_iterator = leafs_iterator + 1;           
                end
                y = [y value];
            end

            %There still may be some branches or leafs left. We will add them
            %below:
            while ( branches_iterator <= length(distance_of_branches) )
                x = [ x distance_of_branches(branches_iterator) ];
                value = value + 1;          
                branches_iterator = branches_iterator + 1;
                y = [y value];
            end
            while ( leafs_iterator <= length(distance_of_leafs) )
                x = [ x , distance_of_leafs(leafs_iterator) ];
                value = value - 1 ;
                leafs_iterator = leafs_iterator + 1;  
                y = [y value];
            end
        end
%%
%%
        function distance_from_rt=distance_from_root(this)
            distance_from_rt = zeros( 1 , this.numNodes );
            root = this.get_root();
                for i=1:this.numNodes
                    node = this.Get_Node(i);
                    dist = distance_of_vectors( node.coordinates , root.coordinates );
                    distance_from_rt(i) = dist;
                end
        end 
%%
%%
%==========================================================================
%> @brief  get_paths_from_leafs_to_root return a vector of vectors of
%vertices. Each individual vector comain vertices in a path from a leaf to
%the root. 
%>
%> @retval vector ov vector of paths between leafs and the root. 
% =========================================================================
    function paths_from_leafs_to_root=get_paths_from_leafs_to_root(this)
        paths_from_leafs_to_root = {};
        %First we need to locate the leafs
        is_leaf_ = this.is_leaf();
        nr_of_paths = 1;
        for i=1:this.numNodes
            if ( is_leaf_(i) )
                %For every leaf we will reconstruct a path from it into the
                %root.
                this_path = [];
                node = this.Get_Node(i);
                while ( node.parent ~= -1 )
                    this_path = [ this_path , node.id ];
                    if ( node.parent ~= -1 ) 
                        node = this.Get_Node( node.parent );
                    end
                end
                %We still want to add the root to the path:
                this_path = [ this_path , node.id ];
                paths_from_leafs_to_root{nr_of_paths} = this_path;
                nr_of_paths = nr_of_paths+1;
            end            
        end
    end 
%%
%%    
%==========================================================================
% This function starts at the leafs and trasverses the tree down to the
% bifrication point and then records the pair ids  of leaf-bifrication with
% all nodes along the path.
%It then take sthe bifrication point and trasverses down to the next 
% bifrication point and so forth untill it hits the soma 
%column1 and colum 2 are the pairs and 
% column 3 --> column n are all the node Id's along the path 

function pairs= node_pairs_with_path_nodes(this)
        is_leaf = this.is_leaf();
        is_branch = this.is_branching_point();
        is_branch(1) =1;
         assignin('base','is_leaf',is_leaf); 
          assignin('base','is_branch',is_branch); 
        nr_of_paths = 1;
        r=1;
         for i=1:this.numNodes

            if ( is_leaf(i) )
               this_path = [];
               node = this.Get_Node(i);
               pairs(r,1)=node.id;
               c=3;
               while ( node.parent ~= -1 )
 
                   if (is_branch(node.parent) == 1)
                     pairs(r,2)=node.parent;                  
                     r=r+1;
                     pairs(r,1)=node.parent;
                     node = this.Get_Node( node.parent );  
                     c=3;
                   else
                     pairs(r,c)=node.parent;
                     node = this.Get_Node( node.parent );
                     
                     c=c+1;
                   end                  
               end 
            end
         end 
        pairs=unique(pairs,'rows');  

        end 
%%
%%
%==========================================================================
% This function starts at the leafs and trasverses the tree down to the
% bifrication point and then records the pair ids  of leaf-bifrication. 
%It then take sthe bifrication point and trasverses down to the next 
% bifrication point and so forth untill it hits the soma 

        function pairs= node_pairs(this)
        is_leaf = this.is_leaf();
        is_branch = this.is_branching_point();
         assignin('base','is_leaf',is_leaf); 
          assignin('base','is_branch',is_branch); 
        nr_of_paths = 1;
        r=1;
         for i=1:this.numNodes
            if ( is_leaf(i) )
               this_path = [];
               node = this.Get_Node(i);
               pairs(r,1)=node.id;
               while ( node.parent ~= -1 )
                   if (is_branch(node.parent) == 1)
                     pairs(r,2)=node.parent;                  
                     r=r+1;
                     pairs(r,1)=node.parent;
                     node = this.Get_Node( node.parent );                     
                   else
                     node = this.Get_Node( node.parent );  
                   end                  
               end 
            end
         end 
        pairs=unique(pairs,'rows');  

        end
%%
%%
    function d=distance_between_two_nodes(node1,node2)
        d=sqrt((node1.coordinates(1)-node2.coordinates(1))^2+(node1.coordinates(2)- node2.coordinates(2))^2+(node1.coordinates(3)-node2.coordinates(3))^2);  
    end 
%%
%%
    function  coordinates = get_coordinates_using_index(this,idx)          
        coordinates=[this.nodesVector(idx).coordinates(1) this.nodesVector(idx).coordinates(2) this.nodesVector(idx).coordinates(3)];
    end 
%%
%%
    function [x_max,y_max, z_max, x_min,y_min,z_min] =get_max_coordinates(this)               
        for i=1:length(this.nodesVector)
          x(i)= this.nodesVector(i).coordinates(1);
          y(i)= this.nodesVector(i).coordinates(2);
          z(i)= this.nodesVector(i).coordinates(3);  
        end 
        x_max=fix(max(x));
        y_max=fix(max(y));
        z_max=fix(max(z));
        x_min=fix(min(x));
        y_min=fix(min(y));
        z_min=fix(min(z));        
    end 
%%
%%    
    function [center_of_mass] =get_center_of_mass(this)
          x=0;
          y=0;
          z=0;
          
        for i=1:length(this.nodesVector)
          x= x +this.nodesVector(i).coordinates(1);
          y=y + this.nodesVector(i).coordinates(2);
          z= z+this.nodesVector(i).coordinates(3);  
        end 
       
        center_of_mass=[x/length(this.nodesVector) y/length(this.nodesVector) z/length(this.nodesVector)]
    end 
%%
%%
%==========================================================================
%> @brief  get_tree_diameter gets the diameter of the tree by calculating
%> the maximum distance from a vector of distances between the soma ad all
%> the nodes in the tree
%>
%> @retval dimeter of the tree
% =========================================================================
    function diameter=get_tree_diameter(this)
       distance_from_rt=this.distance_from_root();
       diameter=max(distance_from_rt);
    end 
%%
%%
    function rad = get_radius(this,rad_type)
               predicate_all = ones(1,length(this.nodesVector));
               predicate_branch = this.is_branching_point();
               predicate_leaf = this.is_leaf();
               predicate=predicate_branch + predicate_leaf; 

               switch rad_type
                    case 'euclidean'
                        euclid_matrix = Spectral_descriptors.Euclidean_distance_matrix(this,predicate_all);
                        rad = zeros(length(predicate),1);
                        for i =1 :length(predicate)
                          if predicate(i)==1
                             rad(i,1) = euclid_matrix(1,i);  
                          end 
                        end 
                        % Remove zero rows
                        rad( all(~rad,2), : ) = [];
                        rad = sort(rad);                    
                    case 'intrinsic'
                        intrinsic_matrix = Spectral_descriptors.intrinsic_distance_matrix( this , predicate_all);   
                        rad = zeros(length(predicate),1);
                        for i =1 :length(predicate)
                          if predicate(i)==1
                             rad(i,1) = intrinsic_matrix(1,i);  
                          end 
                        end 

                        rad( all(~rad,2), : ) = [];
                        rad = sort(rad);    
                    case 'uniform'
                        diameter=this.get_tree_diameter();
                        parts = 30;
                        rad = zeros(length(parts),1);
                        for i =1:parts
                        rad(i,1) = i*(diameter/parts);
                        end  
                        rad = sort(rad);
                    otherwise
                        disp('Not a valid radius type')
               end      
    end 
%%
%%
    function cp=copytree(this)
        cp = copyElement(this);
    end
%%

end  % end public methods

%%
 methods (Access = protected) 
% ======================================================================
%> @brief  THIS FUCNTION CAN BE USED TO COPY ANY HANDLE OBJECTS
%> copyElement method overrides the
%> copyElement@matlab.mixin.Copyable which makes shallow copies. Tree 
%> object stores a handle object (Node) in its properties. If we use simply
%> use copyElement@matlab.mixin.Copyable(obj); to copy the value of Tree that 
%> contains a handle object Node, we are actually copying the handle, not 
%> the object itself. Tree object has properties of type node which are
%> properly cloned via this function. Issue is with the nodesVector which is
%> an array which uuses pointer arithmetic to access array elements.
%> nodeArray object points to the first index of the array and this function
%> is only copying the first index and not the whole array. Tried to remove
%> nodesVector from tree properties and use a function to get nodes using a
%> Queue but performace was drastically affected.
%> This metod is protected. Access from methods in class or subclasses 
%>
%> @param obj the object that needs to be copied
%
%> @retval cp the new object that carries the copy
% ======================================================================

      function cp = copyElement(obj)
%        %copy all elements
         cp = copyElement@matlab.mixin.Copyable(obj);
         % Get handle from root
         hobj = obj.root;         
         % Create default object using the class name 
         new_hobj = eval(class(hobj));         
         % Add public property values from orig object
         Tree.propValues(new_hobj,hobj);
         % Assign the new object to property
         cp.root  = new_hobj;       
      end
end   
%%
%%

 methods (Static)
%==========================================================================
%> @brief  propValues copies the properties between two objects
%>
%> This metod is Static.
%>
%> @param newObj the new copies of the properties
%> @param orgObj the original copies of the properties
% =========================================================================  
      function propValues(newObj,orgObj)
         pl = properties(orgObj);
         for k = 1:length(pl)
            if isprop(newObj,pl{k})
               newObj.(pl{k}) = orgObj.(pl{k});
            end
         end
      end
    
 %%   
end    
end

