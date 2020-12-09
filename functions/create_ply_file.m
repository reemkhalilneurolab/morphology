% ======================================================================
%> @file export_swc_to_ply.m
%> 
%> @brief creates a ply file from matrix of coordinates. ply file is
%>        suitable for vieweing in meshlab
%
%> @param xyz_vector a vector of x,y,z coordinates
%
%> @retval pc ply file object
% ======================================================================

function pc = create_ply_file(xyz_mat)
x=xyz_mat(:,1);
y=xyz_mat(:,2);
z=xyz_mat(:,3);
%> create a point cloud object
pc = pointCloud( [x(:), y(:), z(:)] );
%> save the point cloud object as a .ply file
pcwrite(pc, 'C:\Users\rkhalil\Documents\GitHub\PARS-TMD\Testing\ply_files\file.ply');
%> view it in matlab
pc = pcread('C:\Users\rkhalil\Documents\GitHub\PARS-TMD\Testing\ply_files\file.ply');
%> view it in matlab
pcshow(pc);
end

