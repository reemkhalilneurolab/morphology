classdef Node  < handle
%%
    properties (Access = public)       
        id; % label for parent
        struct_id;
        leftChild=NaN; %left child of the node
        rightSibling=NaN; % right sibling of the left node
        parent;
        coordinates=[];
        radius;
    end
%%
%%
   	methods (Access = public)
%%
%%        
        function new_node = copy_node(this)
            vect_to_copy = [this.id this.struct_id this.coordinates this.radius this.parent ];
            new_node = Node(vect_to_copy);
        end
%%
%%        
 		function this = Node( vect )
            if nargin >0
           %Sample# StructureIdentifier	X Y	Z Radius parentsample     
           this.id = vect(1);
           this.struct_id = vect(2); %1=soma /// 2= axon /// 3=basal dendrite /// 4=apical dendrite 
           this.coordinates = vect(3:5) ;
           this.radius=vect(6);
           this.parent=vect(7);  
            else
            end 
        end 
%%
%%        
         function bool = isnan(this)
            bool = isnan(this.id);
         end
%%
%%
        
    end
end