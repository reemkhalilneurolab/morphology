function ntree= normalize_tree(mytree)
% for taper rate take the max path distance
% pairs= mytree.node_pairs();
% nodesVector = mytree.get_nodesVector();
% 
%    for j=1:size(pairs,1)
%     if pairs(j,1) ~=0 && pairs(j,2)~=0 
%       pairs(j,3)= distance_of_vectors( nodesVector(pairs(j,1)).coordinates , nodesVector(pairs(j,2)).coordinates) ;     
%     end
%     end 


diameter=mytree.get_tree_diameter();
ntree=mytree.copytree();
% ntree=copy_tree(mytree);
num_nodes=mytree.get_num_nodes;
ntree.root.coordinates(1)= mytree.root.coordinates(1)/diameter;
ntree.root.coordinates(2)= mytree.root.coordinates(2)/diameter;
ntree.root.coordinates(3)= mytree.root.coordinates(3)/diameter;
for i=1:num_nodes
  ntree.nodesVector(i).coordinates(1)= mytree.nodesVector(i).coordinates(1)/diameter;
  ntree.nodesVector(i).coordinates(2)= mytree.nodesVector(i).coordinates(2)/diameter;
  ntree.nodesVector(i).coordinates(3)= mytree.nodesVector(i).coordinates(3)/diameter;
end 


% function ntree= normalize_tree(mytree)
% diameter=mytree.get_tree_diameter();
% ntree=mytree.copytree();
% % ntree=copy_tree(mytree);
% num_nodes=mytree.get_num_nodes;
% ntree.root.coordinates(1)= mytree.root.coordinates(1)/diameter;
% ntree.root.coordinates(2)= mytree.root.coordinates(2)/diameter;
% ntree.root.coordinates(3)= mytree.root.coordinates(3)/diameter;
% for i=1:num_nodes
%   ntree.nodesVector(i).coordinates(1)= mytree.nodesVector(i).coordinates(1)/diameter;
%   ntree.nodesVector(i).coordinates(2)= mytree.nodesVector(i).coordinates(2)/diameter;
%   ntree.nodesVector(i).coordinates(3)= mytree.nodesVector(i).coordinates(3)/diameter;
% end 

end

