
function [mytree, all_nodes]=random_walkers(starting_point, probability_to_bifurcate,  probability_to_terminate, number_of_agents, create_simple_move)
    
    mytree=Tree();
    
    all_nodes = [];
    nodes_to_consider = [];
    maximal_number_of_steps = 500;

    
    %Create the root / soma
    %Sample# StructureIdentifier	X Y	Z Radius parentsample
    root = Node([1 1 starting_point(1) starting_point(2) starting_point(3) 1 -1]);
    mytree.addNode( root );
    all_nodes = [ all_nodes , root];
    %Class of an agent:
    
   
    global node_number;
    node_number = 1;
  
    for i = 1:number_of_agents
        % the node has parent -1 and there is no other node with -1 ID . in
        % this case the nodes created here 
%       created a simple move and added a node instead of calling find_next_node_or_nodes as it might create nodes with [0,0,0]
        [new_x,new_y,new_z]=create_simple_move.move_to_next_point(root,mytree);
%         new_node = Node([node_number 1 starting_point(1) starting_point(2) starting_point(3) 1 1]);
        node_number = node_number+1;
        new_node = Node([node_number 1 new_x new_y new_z 1 1]);
        mytree.addNode(new_node);
        
        nodes_to_consider = [ nodes_to_consider , new_node ];
        all_nodes = [ all_nodes , new_node ];
    end

    %This loop will iterate as long as there is alive agent. That will be
    %checked close to the line anyone_there = 0

 
    number_of_steps = 0;
    while ( true )
%         disp( number_of_steps )
        number_of_steps = number_of_steps+1;
        %for every agent:

        new_nodes_to_consider = [];
        for i = 1 : length(nodes_to_consider)      
            v = find_next_node_or_nodes( nodes_to_consider(i)  );
            new_nodes_to_consider = [new_nodes_to_consider,v];
            all_nodes = [ all_nodes , v ];
            number_of_steps = number_of_steps+1;         
        end
        nodes_to_consider = new_nodes_to_consider;

        if ( isempty(nodes_to_consider) )
            break;
        end
        if ( number_of_steps >= maximal_number_of_steps )
            break;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     PrintTree(all_nodes );
  
%      disp(mytree.get_num_nodes());




    function [array_of_nodes]=find_next_node_or_nodes( node )
    %check if it is going to terminate:
        if ( rand() < probability_to_terminate )
            %disp( "terminate" )
            %this agent have died
            array_of_nodes = [];       
        else
            %check if the agent will bifurcate:
            if ( rand() < probability_to_bifurcate )
                %disp( "bifurcate" )
                %the agent will bifurcate. That boils down to add one more
                %agent in the position of the current one. 
                %Here we need to create two new nodes            
                node_number = node_number + 1;
                [new_x,new_y,new_z]=create_simple_move.move_to_next_point(node,mytree);
                node1 = Node([node_number 1 new_x new_y new_z 1 node.id]);
%                 node1 = Node([node_number 1 node.coordinates(1) node.coordinates(2) node.coordinates(3) 1 node.id]);
                mytree.addNode( node1 );
                node_number = node_number + 1;
                [new_x,new_y,new_z]=create_simple_move.move_to_next_point(node,mytree);
%               node2 = Node([node_number 1 node.coordinates(1) node.coordinates(2) node.coordinates(3) 1 node.id]);
                node2 = Node([node_number 1 new_x new_y new_z 1 node.id]);
                mytree.addNode( node2 );
                array_of_nodes = [node1 node2];       
            else
          
                [new_x,new_y,new_z]=create_simple_move.move_to_next_point(node,mytree);
                %the agent will move sample where.
                %disp( "move" )
                %Use this code if no history is to be taken into account:
                node_number = node_number + 1;
                node1 = Node([[node_number 1 new_x new_y new_z 1 node.id]]);
                mytree.addNode( node1 );
                array_of_nodes = [node1];       
            end
        end
    end

assignin('base','all_nodes',all_nodes); 
end

