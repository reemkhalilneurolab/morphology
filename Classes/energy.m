classdef energy
   
    methods(Static)
        
        %==================================================================
        %>@brief   
        %>@param 
        %>@retval
        %==================================================================        
        function get_energy()
            load(fullfile(pwd, CreateTree.save_directory('trees_vector.mat')));
            energy_dist_matrix= energy.get_energy_fun_r_dist_matrix(trees_vector);
            %            energy_dist_matrix= DM.get_energy_dist_matrix(trees_vector);
            assignin('base','energy_dist_matrix',energy_dist_matrix);
            save(fullfile( pwd , CreateTree.save_directory('energy_dist_matrix')),'energy_dist_matrix' );
%             DM.get_clusters(energy_dist_matrix);
            get_dendrogram(energy_dist_matrix,'energy');
        end 
        
        %==================================================================
        %>@brief   
        %>@param 
        %>@retval
        %==================================================================       
        function energy_dist_matrix= get_energy_fun_r_dist_matrix(trees_vector)
            vect_energy =cell(1,length(trees_vector)); 
            vect_energy{1,length(trees_vector)} = [];
            area_energy = zeros(length(trees_vector),2);
%             load(fullfile(pwd, CreateTree.save_directory('vect_energy.mat')));
            for i=1:length(trees_vector)
                disp(['processing Energy for tree :   ,' num2str(i), ' out of ' ,num2str(length(trees_vector)), ' --started at ', datestr(now,'HH:MM:SS.FFF')]);
                tree=normalize_tree(trees_vector(i));
                energy_value= energy.get_energy_func_r(tree);
                vect_energy{i}=energy_value; 
            end 
            save(fullfile( pwd , CreateTree.save_directory('vect_energy')),'vect_energy' );
            energy_dist_matrix=zeros(length(trees_vector),length(trees_vector));
            for i=1:length(trees_vector)
                 for j=(i+1):length(trees_vector)
                       denergy=DM.get_distance_of_plots_grid_version( vect_energy{i}(:,1),vect_energy{i}(:,2),vect_energy{j}(:,1),vect_energy{j}(:,2) );
                       energy_dist_matrix(i,j) = denergy;
                       energy_dist_matrix(j,i) = denergy;
                 end 
                area_energy(i,1) = trapz(vect_energy{i}(:,1),vect_energy{i}(:,2));
                if vect_energy{i}(1,1) == 0 && vect_energy{i}(1,1) ~= 0
                    area_energy(i,2) = vect_energy{i}(1,2);
                else
                    area_energy(i,2) = vect_energy{i}(end,2);               
                end 
                 
            end  
            save(fullfile( pwd , CreateTree.save_directory('area_energy')),'area_energy' );
 
        end 
        
        %==================================================================
        %>@brief function gets the energy  contribution of every node to the
        %> soma. 
        %>@param 
        %>@retval
        %==================================================================        
        function energy_value= get_energy_func_r(mytree)
% 
            is_leaf = mytree.is_leaf( );
            is_branch = mytree.is_branching_point( );
            
            positions_and_diam = [];
            root = mytree.get_root();
            pos = [root.coordinates(1) root.coordinates(2) root.coordinates(3)];
            for i=1:mytree.numNodes
               
               if ( (is_leaf(i) == 1) || (is_branch(i) == 1) )
                   node = mytree.Get_Node(i);
                   v1 =  [node.coordinates(1) node.coordinates(2) node.coordinates(3)];
                   v2= distance_of_vectors( v1 , pos );
                   v = [v2 v1];
%                    v2 = [ v node.radius ]; 
                   positions_and_diam = [ positions_and_diam ; v ];
               end
            end 
        
            
            energy_value = 0;
            contribution_1=0;
            contribution_2=0;
            contribution_3=0;
            positions_and_diam = sortrows(positions_and_diam,1);
     
            for i = 1:size(positions_and_diam)        
                %I am adding +1 to denominator to make sure that we are not
                %dividing by zero if by any chance we get a point very
                %close to the root. 
%                 dist = distance_of_vectors( positions_and_diam(i,1:3),pos );
                dist=positions_and_diam(i,1);
%                 contribution_term_1 = positions_and_diam(i,4) * (pos(1)-positions_and_diam(i,1)) / ( sqrt(dist + 1) );
%                 contribution_term_2 = positions_and_diam(i,4) * (pos(2)-positions_and_diam(i,2)) / ( sqrt(dist + 1) );
%                 contribution_term_3 = positions_and_diam(i,4) * (pos(3)-positions_and_diam(i,3)) / ( sqrt(dist + 1) );
%                 contribution_term_1 = (1/positions_and_diam(i,4)) * (pos(1)-positions_and_diam(i,1)) / ( sqrt(dist) );
%                 contribution_term_2 = (1/positions_and_diam(i,4)) * (pos(2)-positions_and_diam(i,2)) / ( sqrt(dist) );
%                 contribution_term_3 = (1/positions_and_diam(i,4)) * (pos(3)-positions_and_diam(i,3)) / ( sqrt(dist) );
%                 
                contribution_term_1 =  (pos(1)-positions_and_diam(i,2)) / ( dist );
                contribution_term_2 =  (pos(2)-positions_and_diam(i,3)) / ( dist );
                contribution_term_3 =  (pos(3)-positions_and_diam(i,4)) / ( dist );
            
               contribution_1=contribution_1+contribution_term_1;
               contribution_2=contribution_2+contribution_term_2;
               contribution_3=contribution_3+contribution_term_3; 
              
               energy_value(i,1)=dist;
               energy_value(i,2)=contribution_1^2+contribution_2^2+contribution_3^2;
               
            end                          
        end
        %==================================================================
        
%         %==================================================================
%         %>@brief function gets the energy  contribution of every node to the
%         %> soma. 
%         %>@param 
%         %>@retval
%         %==================================================================        
%         function energy_value= get_energy_func_r(mytree)
% % 
%             is_leaf = mytree.is_leaf( );
%             is_branch = mytree.is_branching_point( );
%             
%             positions_and_diam = [];
%             for i=1:mytree.numNodes
%                
%                if ( (is_leaf(i) == 1) || (is_branch(i) == 1) )
%                    node = mytree.Get_Node(i);
%                    v = [ node.coordinates(1) node.coordinates(2) node.coordinates(3) node.radius ]; 
%                    positions_and_diam = [ positions_and_diam ; v ];
%                end
%             end 
%         
%             
%             energy_value = 0;
%             sizes = size(positions_and_diam);
%             contribution_1=0;
%             contribution_2=0;
%             contribution_3=0;
%             root = mytree.get_root();
%              pos = [root.coordinates(1) root.coordinates(2) root.coordinates(3)];
% %             pos =  mytree.get_center_of_mass();
%             dist1 = zeros(size(positions_and_diam));
%             for i=1:size(positions_and_diam)
%             dist1(i) = distance_of_vectors( positions_and_diam(i,1:3),pos );
%             end
%             dist1 = sortrows(dist1,1);
%      
%             for i = 1:sizes(1)        
%                 %I am adding +1 to denominator to make sure that we are not
%                 %dividing by zero if by any chance we get a point very
%                 %close to the root. 
% %                 dist = distance_of_vectors( positions_and_diam(i,1:3),pos );
%                 dist=dist1(i);
% %                 contribution_term_1 = positions_and_diam(i,4) * (pos(1)-positions_and_diam(i,1)) / ( sqrt(dist + 1) );
% %                 contribution_term_2 = positions_and_diam(i,4) * (pos(2)-positions_and_diam(i,2)) / ( sqrt(dist + 1) );
% %                 contribution_term_3 = positions_and_diam(i,4) * (pos(3)-positions_and_diam(i,3)) / ( sqrt(dist + 1) );
% %                 contribution_term_1 = (1/positions_and_diam(i,4)) * (pos(1)-positions_and_diam(i,1)) / ( sqrt(dist) );
% %                 contribution_term_2 = (1/positions_and_diam(i,4)) * (pos(2)-positions_and_diam(i,2)) / ( sqrt(dist) );
% %                 contribution_term_3 = (1/positions_and_diam(i,4)) * (pos(3)-positions_and_diam(i,3)) / ( sqrt(dist) );
% %                 
%                 contribution_term_1 =  (pos(1)-positions_and_diam(i,1)) / ( dist );
%                 contribution_term_2 =  (pos(2)-positions_and_diam(i,2)) / ( dist );
%                 contribution_term_3 =  (pos(3)-positions_and_diam(i,3)) / ( dist );
%             
%                contribution_1=contribution_1+contribution_term_1;
%                contribution_2=contribution_2+contribution_term_2;
%                contribution_3=contribution_3+contribution_term_3; 
%               
%                energy_value(i,1)=dist;
%                energy_value(i,2)=contribution_1^2+contribution_2^2+contribution_3^2;
%                
%             end                          
%         end
%         %==================================================================
        %>@brief    
        %In this function we compute the value of a energy at the point pos.
        %The input is the collection of nodes of the tree and the diameters
        %of the tree therein stored in positions_and_diam. 
        %>@param positions_and_diam - (x,y,z,r), where x,y,z are coordinates
        %of nodes and r is the radius of the tree threin
        %>@param pos - vector of three elements describing (x,y,z) position
        %of a point. 
        %>@retval Value of the energy of a tree at a given point.
        %================================================================== 
        function energy_value = compute_energy_at_given_point( positions_and_diam , pos )
            energy_value = 0;
            sizes = size(positions_and_diam);
            contribution_1=0;
            contribution_2=0;
            contribution_3=0;
            for i = 1:sizes(1)        
                contribution_term_1 = (positions_and_diam(i,4) * pos(1)-positions_and_diam(i,1)) / (distance_of_vectors( positions_and_diam(i,1:3),pos )  );
                contribution_term_2 = (positions_and_diam(i,4) * pos(2)-positions_and_diam(i,2)) / (distance_of_vectors( positions_and_diam(i,1:3),pos )  );
                contribution_term_3 = (positions_and_diam(i,4) * pos(3)-positions_and_diam(i,3)) / (distance_of_vectors( positions_and_diam(i,1:3),pos )  );

               contribution_1=contribution_1+contribution_term_1;
               contribution_2=contribution_2+contribution_term_2;
               contribution_3=contribution_3+contribution_term_3;        
            end

            energy_value=contribution_1^2+contribution_2^2+contribution_3^2;
        end
        
        %==================================================================
        %>@brief    
        %In this function we compute the values of a energy function at a
        %collection of 3d points stored in the vector points_vector.
        %>@param tree , vector of points in R^3 (having three coordinates) e.g. rand(10,3). 
        %>@retval Vector of valuse of energy of a mytree tree computed at the
        %points in points_vector.
        %==================================================================
        function energy_for_points = compute_energy_for_points( mytree , points_vector )
            is_leaf = mytree.is_leaf( );
            is_branch = mytree.is_branching_point( );
            positions_and_diam = [];
            for i=1:mytree.numNodes
               node = mytree.Get_Node(i);
               %if it is leaf or brach
               if ( (is_leaf(i) == 1) || (is_branch(i) == 1) )
                   %We get its coodrginates:
                   v = [ node.coordinates(1) node.coordinates(2) node.coordinates(3) node.radius ]; 
                   positions_and_diam = [ positions_and_diam ; v ];               
               end
            end

            energy_for_points = zeros(1,length(points_vector));
            %For every point in the points_vector, compute the value of its energy
            %given positions_and_diam vector. 
            for pt = 1:length(points_vector)
                pos = points_vector(pt,:);
                 energy_for_points(pt) = energy.compute_energy_at_given_point( positions_and_diam , pos );                       
            end
        end

        %==================================================================
        %>@brief THIS FUNCTION IS NOT USED ..COMPUTATIONALY EXPENSIVE   
        %In this function we compute the value of a energy at grid of points.
        %>@param tree , min and max x,y and z values of the grid as well as
        %dx which is a diameter of a tree. 
        %>@retval Value of the energy of a tree at a grid of points. 
        % this fucntion will not work with the
        % compute_energy_for_pointssince we do not have a grid and hence a
        % volume. 
        %==================================================================
%         function energy_on_grid = compute_energy_on_grid( mytree , min_x , max_x , min_y , max_y , min_z , max_z , dx )  
%                 is_leaf = mytree.is_leaf( );
%                 is_branch = mytree.is_branching_point( );
% 
%                 positions_and_diam = [];
%                 for i=1:mytree.numNodes
%                    node = mytree.Get_Node(i);
%                    if ( (is_leaf(i) == 1) || (is_branch(i) == 1) )
%                        v = [ node.coordinates(1) node.coordinates(2) node.coordinates(3) node.radius ]; 
%                        positions_and_diam = [ positions_and_diam ; v ];
%                    end
%                 end
% 
%                 number_of_iterations_x = (max_x-min_x)/dx;
%                 number_of_iterations_y = (max_y-min_y)/dx;
%                 number_of_iterations_z = (max_z-min_z)/dx;
% 
%                 energy_on_grid = zeros( number_of_iterations_x , number_of_iterations_y , number_of_iterations_z );
%                 for i = 1:number_of_iterations_x
%                     for j = 1:number_of_iterations_y
%                         for k = 1:number_of_iterations_z
%                             pos = [ min_x+i*dx , min_y + j*dx , min_z+k*dx ];
%                             energy_on_grid(i,j,k) = energy.compute_energy_at_given_point( positions_and_diam , pos );                
%                         end
%                     end
%                 end
%         end

        %==================================================================
        %>@brief    
        %In this function we compute the value of a energy at grid of points.
        %>@param tree , min and max x,y and z values of the grid as well as
        %dx which is a diameter of a tree. 
        %>@retval Value of the energy of a tree at a grid of points. 
        % this fucntion will not work with the
        % compute_energy_for_pointssince we do not have a grid and hence a
        % volume. 
        %==================================================================
        function energy_on_grid = compute_energy_on_grid( mytree , min_x , max_x , min_y , max_y , min_z , max_z , dx )  
            is_leaf = mytree.is_leaf( );
            is_branch = mytree.is_branching_point( );

            positions_and_diam = [];
            for i=1:mytree.numNodes
               node = mytree.Get_Node(i);
               if ( (is_leaf(i) == 1) || (is_branch(i) == 1) )
                   v = [ node.coordinates(1) node.coordinates(2) node.coordinates(3) node.radius ]; 
                   positions_and_diam = [ positions_and_diam ; v ];
               end
            end

%             number_of_iterations_x = abs((max_x-abs(min_x))/dx);
%             number_of_iterations_y = abs((max_y-abs(min_y))/dx);
%             number_of_iterations_z = abs((max_z-abs(min_z))/dx);
            number_of_iterations_x = 100;
            number_of_iterations_y = 100;
            number_of_iterations_z = 100;
            
            N=number_of_iterations_x*number_of_iterations_y*number_of_iterations_z;
            K=150;
            c = randsample(N,K);
            
            energy_on_grid = zeros( number_of_iterations_x , number_of_iterations_y , number_of_iterations_z );
            for i = 1:number_of_iterations_x
                for j = 1:number_of_iterations_y
                    for k = 1:number_of_iterations_z
                       res = ismember(i*j*k,c);
                       if res==1
                        pos = [ min_x+i*dx , min_y + j*dx , min_z+k*dx ];
                        energy_on_grid(i,j,k) = energy.compute_energy_at_given_point( positions_and_diam , pos );
                        
                       else
                         energy_on_grid(i,j,k)=0;  
                       end
                       
                    end
                end
            end
         
        end

        %==================================================================
        %>@brief   
        %>@param 
        %>@retval
        %==================================================================
%         function test_energy_volume()
%             load('C:\Users\rkhalil\Documents\GitHub\PARS-TMD\test_cases\energy_on_grid_.mat');
%             load('C:\Users\rkhalil\Documents\GitHub\PARS-TMD\test_cases\trees_vector\trees_vector.mat');
%             for i=1:length(trees_vector)
%             min_energy= 0;
%             max_energy =1000000;
%             vect_energy = energy.energetic_volume( energy_on_grid_{i} , min_energy, max_energy )
%             energetic_volume(i)=vect_energy;
%             end 
%             save(fullfile( pwd , "\test_cases\energetic_volume" ),'energetic_volume' );
%         end 
        
        %==================================================================
        %>@brief    
        %This function computes the plot of the values of energy of points 
        %in the energy_on_grid in between min_energy and max_energy.
        %>@param energy on a grid, min and max value of energy to be taken
        %into consideration
        %>@retval A plot of energy values of point having energy values in
        %between min_energy and max_energy. There length of the returned
        %array is equal to the number of points having the intermediate
        %energy levels. 
        %==================================================================
        function energetic_volume_ = energetic_volume_plot( energy_on_grid , min_energy, max_energy )
            nr = numel( energy_on_grid );
            vect = reshape( energy_on_grid , [1,nr] );
            vect = sort( vect );
            %now we need to extract part between min_energy and max_energy.
            %To do this, we will find positions in which the minimal and
            %maximal energy occures in the vector. The beginning posiiton
            %will be indicated by it1, the end by it2. 
            it1 = 1;
            while ( (it1 < length(vect) ) && (vect(it1) < min_energy) )
                it1 = it1 + 1;
            end
            it2 = it1;
            while ( (it2 < length(vect)) && (vect(it2) < max_energy) )
                it2 = it2+1;
            end
            
            energetic_volume_ = [];
            if ( it1 ~= length(vect) )
            
                energetic_volume_ = zeros( 2,it2-it1 );
                energetic_volume_(1,1) = it1;
                energetic_volume_(1,2) = vect(it1);

                for i=(it1+1):it2
                    energetic_volume_(1,i-it1) = vect(i);
                    energetic_volume_(2,i-it1) = i/nr;
                end
            end        
        end
        
        %==================================================================
        %>@brief    
        %This function computes the percentage of a volume of the whole
        %grid occupied by points with enedgy levels beteeen min_energy and
        %max_energy.
        %>@param energy on a grid, min and max value of energy to be taken
        %into consideration
        %>@retval Percentage of the volume occupuef by points of the grid
        %ahving energy valyes in between min_energy and max_energy.
        %==================================================================
        function energetic_volume = energetic_volume( energy_on_grid , min_energy, max_energy )
            energetic_volume_plt = energy.energetic_volume_plot( energy_on_grid , min_energy, max_energy );
            energetic_volume = length( energetic_volume_plt )/numel( energy_on_grid );
        end      

        %==================================================================
        %>@brief  NOT USED ..EXPENSIVE 
        %>@param 
        %>@retval
        %==================================================================  
        function energy_dist_matrix= get_energy_dist_matrix1(trees_vector)
          for i=1:length(trees_vector)
               mytree=trees_vector(i);
               [x_max,y_max, z_max, x_min,y_min,z_min] =mytree.get_max_coordinates();
               energy_on_grid = energy.compute_energy_on_grid( mytree , x_min , x_max , y_min , y_max , z_min , z_max , 0.1 );
               %Here is the plot of enery levels in between extrema. If we want anything
               %inside the range, those numbers have to be chosen carefully, as i did
               %here. If not, just put 0, very large number and that should do.
               %this is to see a part:
               %ew = energy.energetic_volume( energy_on_grid , 2000, 3000 );
               %this is to see the whole thing:
               energy_vector = energy.energetic_volume_plot( energy_on_grid , 0, 1000000 );
               assignin('base','energy_on_grid',energy_on_grid);
               assignin('base','energy_vector',energy_vector);
               vect_energy{i}=energy_vector;
               assignin('base','vect_energy',vect_energy); 
        %                plot(ew(1,:),ew(2,:));
               energy_on_grid_{i}=energy_on_grid;

          end 
          save(fullfile( pwd , CreateTree.save_directory('energy_on_grid_')),'energy_on_grid_' );
          for i=1:length(trees_vector)
                 for j=(i+1):length(trees_vector)
                       denergy=DM.get_distance_of_plots_grid_version( vect_energy{i}(:,1),vect_energy{i}(:,2),vect_energy{j}(:,1),vect_energy{j}(:,2) ); 
                       energy_dist_matrix(i,j) = denergy;
                       energy_dist_matrix(j,i) = denergy;
                 end 
          end
        end 

    end
end

