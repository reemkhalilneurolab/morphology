classdef GRAPH
    
   methods(Static)
       
        function G = convert_to_matlab_graph(mytree) 
            nodes = Queue('Node') ;      
            nodes.add( mytree.get_root() );
            s = [];
            t = [];
            weights = [];
            while ( ~nodes.isempty() )
                node= nodes.pop();
                children = mytree.get_children(node);
                for i = 1:length(children)
                    s = [s,node.id];
                    t = [t,children(i).id];
                    weights = [weights distance_of_vectors( node.coordinates , children(i).coordinates )];
                    nodes.add( children(i) );
                end              
            end
            G = graph(s,t,weights);
        end  
        
        function G = convert_to_matlab_graph_all_edges_same_distance(mytree) 
            nodes = Queue('Node') ;      
            nodes.add( mytree.get_root() );
            s = [];
            t = [];
            while ( ~nodes.isempty() )
                node= nodes.pop();
                children = mytree.get_children(node);
                for i = 1:length(children)
                    s = [s,node.id];
                    t = [t,children(i).id];
                    nodes.add( children(i) );
                end              
            end
            G = graph(s,t);
        end  
  
    end
end

