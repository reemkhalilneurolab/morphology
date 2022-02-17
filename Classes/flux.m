classdef flux
 

    methods(Static)
        
        %==================================================================
        %>@brief   This descriptor associates to a given radius r the sum 
        %>         of the angles between dendrites and normal directions to the sphere
        %>         of radius r centered at the soma, at the points where the dendrites
        %>         intersect with that sphere.
        %>         we  subtract the number of leaves, all within a given radius r from the soma. 
        %>@param none
        %>@retval generates a dendrogram and a distance matrix
        %==================================================================        
         function  get_flux() 
            load(fullfile(pwd, CreateTree.save_directory('trees_vector.mat')));
    
            [flux_dist_matrix, north_flux_dist_matrix, south_flux_dist_matrix, lateral_flux_dist_matrix] = flux.get_dist_matrix(trees_vector);
        
            assignin('base','flux_dist_matrix',flux_dist_matrix);
            assignin('base','north_flux_dist_matrix',north_flux_dist_matrix);
            assignin('base','south_flux_dist_matrix',south_flux_dist_matrix);
            assignin('base','lateral_flux_dist_matrix',lateral_flux_dist_matrix);
            
            save(fullfile( pwd , CreateTree.save_directory('flux_dist_matrix')),'flux_dist_matrix' );
            save(fullfile( pwd , CreateTree.save_directory('north_flux_dist_matrix')),'north_flux_dist_matrix' );
            save(fullfile( pwd , CreateTree.save_directory('south_flux_dist_matrix')),'south_flux_dist_matrix' );
            save(fullfile( pwd , CreateTree.save_directory('lateral_flux_dist_matrix')),'lateral_flux_dist_matrix' );

            get_dendrogram(flux_dist_matrix,'flux');
            get_dendrogram(north_flux_dist_matrix,'North_Flux');
            get_dendrogram(south_flux_dist_matrix,'South_Flux');
            get_dendrogram(lateral_flux_dist_matrix,'Lateral_Flux');

        end

        function [flux_dist_matrix ,north_flux_dist_matrix, south_flux_dist_matrix, lateral_flux_dist_matrix] = get_dist_matrix(trees_vector)
            area_flux = zeros(length(trees_vector),2);
            vect_flux = cell(1,length(trees_vector));
            vect_flux{1,length(trees_vector)} = [];  
            
            vect_south = cell(1,length(trees_vector));
            vect_south{1,length(trees_vector)} = [];
            
            vect_north = cell(1,length(trees_vector));
            vect_north{1,length(trees_vector)} = [];
            
            vect_Later = cell(1,length(trees_vector));
            vect_Later{1,length(trees_vector)} = [];
            
            
            for i=1:length(trees_vector)
                 tree=normalize_tree(trees_vector(i));
                 disp(['processing flux for tree :   ,' num2str(i), ' out of ' ,num2str(length(trees_vector)), ' --started at ', datestr(now,'HH:MM:SS.FFF')]);
                
                [f ,n, s, L]=flux.get_flux_sholl_descriptor(tree); 
                vect_north{i} = n;
                vect_south{i} = s;
                vect_Later{i} = L;
                vect_flux{i} =  f;
                
            end 
           save(fullfile( pwd , CreateTree.save_directory('vect_flux')),'vect_flux' ); 
           save(fullfile( pwd , CreateTree.save_directory('vect_north')),'vect_north' );
           save(fullfile( pwd , CreateTree.save_directory('vect_south')),'vect_south' );
           save(fullfile( pwd , CreateTree.save_directory('vect_Later')),'vect_Later' );           
            
            flux_dist_matrix = zeros(length(trees_vector),length(trees_vector));
            north_flux_dist_matrix = zeros(length(trees_vector),length(trees_vector));
            south_flux_dist_matrix = zeros(length(trees_vector),length(trees_vector));
            lateral_flux_dist_matrix = zeros(length(trees_vector),length(trees_vector));           
     
            
            for i=1:length(trees_vector)
                 for j=(i+1):length(trees_vector)   
                   dflux=DM.get_distance_of_plots_grid_version( vect_flux{i}(:,1),vect_flux{i}(:,2),vect_flux{j}(:,1),vect_flux{j}(:,2) ); 
                   dnorth=DM.get_distance_of_plots_grid_version( vect_north{i}(:,1),vect_north{i}(:,2),vect_north{j}(:,1),vect_north{j}(:,2) );
                   dsouth=DM.get_distance_of_plots_grid_version( vect_south{i}(:,1),vect_south{i}(:,2),vect_south{j}(:,1),vect_south{j}(:,2) );
                   dlateral=DM.get_distance_of_plots_grid_version( vect_Later{i}(:,1),vect_Later{i}(:,2),vect_Later{j}(:,1),vect_Later{j}(:,2) );
            
                   flux_dist_matrix(i,j) = dflux;
                   flux_dist_matrix(j,i) = dflux;
                   north_flux_dist_matrix(i,j) = dnorth;
                   north_flux_dist_matrix(j,i) = dnorth;
                   south_flux_dist_matrix(i,j) = dsouth;
                   south_flux_dist_matrix(j,i) = dsouth;
                   lateral_flux_dist_matrix(i,j) = dlateral;
                   lateral_flux_dist_matrix(j,i) = dlateral;
                 end
                area_flux(i,1) = trapz(vect_flux{i}(:,1),vect_flux{i}(:,2));
                if vect_flux{i}(1,1) == 0 && vect_flux{i}(1,2) ~= 0
                    area_flux(i,2) = vect_flux{i}(1,2);
                else
                    area_flux(i,2) = vect_flux{i}(end,2);               
                end 
                 
                [m , b ] = get_linear_regression(vect_flux{i}(:,1),vect_flux{i}(:,2));
                area_flux(i,3) = m;
                area_flux(i,4) = b; 
                 
            end
            save(fullfile( pwd , CreateTree.save_directory('area_flux')),'area_flux' );

        end 

        function [flux1, northFlxSum, southFlxSum, lateralFluxSum] = get_flux_sholl_descriptor(mytree)
            distance_from_root = zeros( 1 , mytree.numNodes );
            root = mytree.get_root();
            %For every branch, find the distance of the midpoint of that branch to the root. 
            rad = [];%           
            flux_mat=[];
            flux1=[];
            flux_mat = flux.calculate_flux(mytree);            
            flux1= [flux1;flux_mat(:,1:2)];
            [somaVect min_y max_y] = flux.get_contour_soma_points(mytree);
            %In the code below the whole flux is divided into north,
            %south and lateral.             
            lateralFlux=flux_mat(flux_mat(:,4)>=min_y & flux_mat(:,4)<=max_y,:);
            % y values determine direction
            northFlux=flux_mat(flux_mat(:,4)>max_y,:);  
            southFlux=flux_mat(flux_mat(:,4)<min_y,:);            
            
            [flx,~,indx]=unique(flux1(:,1)); 
            flux1 = [flx,accumarray(indx,flux1(:,2),[],@sum)];

            %North Flux---add fluxes at different radii
            [flx1,~,indx1]=unique(northFlux(:,1)); 
            northFlxSum = [flx1,accumarray(indx1,northFlux(:,2),[],@sum)];

            %South Flux---add fluxes at different radii
            [flx2,~,indx2]=unique(southFlux(:,1)); 
            southFlxSum = [flx2,accumarray(indx2,southFlux(:,2),[],@sum)];

            %Lateral Flux between min_y and max_y
            [flx3,~,indx3]=unique(lateralFlux(:,1)); 
            lateralFluxSum = [flx3,accumarray(indx3,lateralFlux(:,2),[],@sum)];
        end 
       
        %==================================================================
        %>@brief   
        %>@param 
        %>@retval
        %==================================================================
        function fluxMatrix = calculate_flux( tree )
            fluxMatrix=[];  
            pwp=tree.node_pairs_with_path_nodes();
            nodesVector = tree.get_nodesVector();
            pridicate_all = ones(1,length(tree.nodesVector));
            predicate_branch = tree.is_branching_point();
            predicate_leaf = tree.is_leaf();
            predicate=predicate_branch + predicate_leaf;

            % get all the radii i.e the distances between the soma to all
            % bifricating points
            rad=[];
            root_coord = tree.get_root().coordinates;
                for i =1 :length(predicate)
                  if predicate(i)==1
                     rad = [rad; distance_of_vectors( nodesVector(i).coordinates , root_coord )];  
                  end 
                end  
                rad= sort(rad);  
                new_rad = [];
                for t=1:size(rad,1)-1
                    new_rad(t) = (rad(t,1)+rad(t+1,1))/2;                    
                end 
                rad = new_rad';
                rad( all(~rad,2), : ) = [];
            %get euclidean distance matrix between all points/nodes/markers
            euclid_matrix = Spectral_descriptors.Euclidean_distance_matrix(tree,pridicate_all);          
            %swap the first two colunns in the path matrix so that
            %when walking along the path we can reach the last point
            %in column 2
            v = [1 2];
            pwp(:, v) = pwp(:, flip(v));
            for i=1:size(rad,1) %loop through radii
                pts_mtx =[];
                for p=1:size(pwp,1) % loop through all the paths
                    %for a given path walk along the path and get the
                    %two points ..one outside and one inside the sphere
                       if (pwp(p,1) && pwp(p,2)~=0) % non of the two points are zeros 

                            submtx = pwp(p,:); %extract the row 
                            submtx( :, all(~submtx,1) ) = [];% delete zero columns
                            if size(submtx,2) > 2 % make sure we have a path
                                 for f=size(submtx,2)-1:-1:2 % walk along the path 
                                    if euclid_matrix(1,submtx(1,f)) >= rad(i) && euclid_matrix(1,submtx(1,f+1)) <= rad(i)
                                       pts_mtx = [pts_mtx ; submtx(1,f) submtx(1,f+1)]; 
                                       break;
                                    end 

                                 end 
                            end                      

                       end 
                end 
                % start processing the flux at raduis i
                for h=1:size(pts_mtx,1)
                  A = nodesVector(pts_mtx(h,1)).coordinates; %child
                  B = nodesVector(pts_mtx(h,2)).coordinates ; %parent
                  fluxVect= [rad(i) flux.linearize(B(1,1),B(1,2),B(1,3),A(1,1),A(1,2),A(1,3),rad(i))];
                  fluxMatrix=[fluxMatrix; fluxVect];   
                end 

                
             end 

        end

        %==================================================================
        %>@brief   The bary center center coordinate of a point of
        %>         intersection between A and B
        %>@param 
        %>@retval
        %==================================================================
        function fluxVect=linearize(x1,y1,z1,x2,y2,z2,rad)
%           x1,y1,z1 is the father
            flux=0;
            %============== point of  Intersection on Sphere ===
            A=(x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2; 
            B=2*((x1*x2)+(y1*y2)+(z1*z2)-x1^2-y1^2-z1^2);
            C= x1^2 + y1^2 + z1^2;

            bSquare=B*B;
            minusB=-1*B;
            fourTimesA=4*A;
            CminusRsquare=C-rad^2;
            twoTimesA=2*A;
            squareroot= sqrt(bSquare- (fourTimesA*CminusRsquare));

            tminus=(minusB-squareroot)/twoTimesA;
            tplus=(minusB+squareroot)/twoTimesA;

            t=tplus;
            %point of intersection
            Cx=(1-t)*x1 + t*x2;
            Cy=(1-t)*y1 + t*y2;
            Cz=(1-t)*z1 + t*z2;
            %==============  FLUX   ==============================
            Nc=1/(sqrt(Cx^2+Cy^2+Cz^2) * sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2 ) ); %the normal vector
            Fc=Cx*(x2-x1) + Cy*(y2-y1) + Cz*(z2-z1) ;
            flux=Nc*Fc;
            fluxVect=[flux Cx Cy Cz];
 
        end
        
        %==================================================================
        %>@brief In this function we obtain the coordinates of the somma. 
        %>       It is used to get the y (or in math language, z) 
        %>       coordinates span of the some. They are used to 
        %>       determine/define the lateral branches in oppose to north 
        %>       and south branches.   
        %>@param 
        %>@retval
        %==================================================================
        function [somaVect min_y max_y] = get_contour_soma_points(mytree)
        %single soma contour that is averaged into N soma points
        somaVect=[]; %to be used later
        %min_y and max_y are the epsilon range values
        min_y=0;
        max_y=0;
         for i = 1:mytree.get_num_nodes()
            node = mytree.Get_Node(i);
            if (node.struct_id==1) %struct_id =1 in the swc file
                somaVect=[somaVect node];
                if node.coordinates(2)< min_y
                    min_y=node.coordinates(2);
                end 
                 if node.coordinates(2)> max_y
                    max_y=node.coordinates(2);
                end
            end
         end
        end %get_contour_soma_points

    end
end

