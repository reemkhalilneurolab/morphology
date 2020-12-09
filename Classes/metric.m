classdef metric 
   
methods(Static)
    function metriclearn()
        tStart = tic; 
        %In this file we will read distance matrices comong from various descriptor
        %and build a combinations of them to get one distance matrix to rule them
        %all.
        load(fullfile(pwd, CreateTree.save_directory('swc_vector.mat')));
            %First we read the distance matrices: 
            load(fullfile(pwd, CreateTree.save_directory('branching_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('tortuosity_dist_matrix.mat')));                
            load(fullfile(pwd, CreateTree.save_directory('flux_dist_matrix.mat')));               
            load(fullfile(pwd, CreateTree.save_directory('leaf_index_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('energy_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('wiring_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('tmd_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('taper_rate_dist_matrix.mat')));
      
        tortuosity      = tortuosity_dist_matrix;        
        branching       = branching_dist_matrix;        
        flux            = flux_dist_matrix;      
        leaf_index      = leaf_index_dist_matrix;        
        energy          = energy_dist_matrix;
        wiring          = wiring_dist_matrix;
        tmd      = tmd_dist_matrix;
         taper      = taper_rate_dist_matrix;
        
      
        %Now we create weights with which they should be used. 
        %original
        dist_mat_1 = {energy flux leaf_index branching wiring tortuosity  tmd taper};

        [r,c] = size(dist_mat_1);
        dist_mat =[];
        j=1;
        for i=1:c
           hasnan = any(isnan(dist_mat_1{i}(:)));
           if hasnan == 0
             dist_mat{j} = dist_mat_1{i};
             j=j+1;
           end
        end
             
        
        [classes,~,indx]=unique(swc_vector(4,:)');        
        num_of_descriptors = length (dist_mat);
        num_of_classes = length (classes);
        [~ , num_of_neurons] = size(swc_vector);
        avg_dist_mtx = [];
        indices_mat=[];
        flattened_mtx =[];
        for j=1:num_of_descriptors 
          mat = cell2mat(dist_mat(j)); 
          avg_distance_within_i_class =[];          
          for i=1:length(classes)
                    % get all the indices of the class i 
                    k = find(indx==i);
                    %starting index
                    s = min(k);
                    %ending index
                    e = max(k); 
%                     indices_mat = [ indices_mat; [i s e]];
                    avg_distance_within_i_class = [avg_distance_within_i_class  mean(mat(1:num_of_neurons,s:e),2)];                                   
                    ahmad = mat(1:num_of_neurons,s:e);         
          end
          
          avg_dist_mtx = [avg_dist_mtx  avg_distance_within_i_class];
%           avg1= avg_dist_mtx.';
          
        end 
        
        for i=1:length(classes)
                    % get all the indices of the class i 
                    k = find(indx==i);
                    %starting index
                    s = min(k);
                    %ending index
                    e = max(k); 
                    indices_mat = [ indices_mat; [i s e]];           
        end
      
      indx1=indx-1;
      save(fullfile( pwd , CreateTree.save_directory('avg_dist_mtx')),'avg_dist_mtx' ); 
      xlswrite(fullfile(pwd,'\save\indices_mat.xlsx'), indices_mat, 'Sheet1', 'A1'); 
      xlswrite(fullfile(pwd,'\save\vectorized_mat_supervised.xlsx'), avg_dist_mtx, 'Sheet1', 'A1');  
      xlswrite(fullfile(pwd,'\save\labels.xlsx'), indx1, 'Sheet1', 'A1'); 
      
                
    end
    
    function vectorize_unsupervised()
        load(fullfile(pwd, CreateTree.save_directory('swc_vector.mat')));
        load(fullfile(pwd, CreateTree.save_directory('area_branching_pattern.mat')));
        load(fullfile(pwd, CreateTree.save_directory('area_tortuosity.mat')));                
        load(fullfile(pwd, CreateTree.save_directory('area_flux.mat')));               
        load(fullfile(pwd, CreateTree.save_directory('area_leaf_index.mat')));
        load(fullfile(pwd, CreateTree.save_directory('area_energy.mat')));
        load(fullfile(pwd, CreateTree.save_directory('area_wiring.mat'))); 
        load(fullfile(pwd, CreateTree.save_directory('area_taper_rate.mat')));

        area_mat = {area_energy area_flux  area_leaf_index  area_branching_pattern  area_wiring  area_tortuosity area_taper_rate };
        vectorized_unsupervised = [];
        for j=1:length (area_mat)
          mat = cell2mat(area_mat(j));
          vectorized_unsupervised = [vectorized_unsupervised mat];          
        end 
      save(fullfile( pwd , CreateTree.save_directory('vectorized_unsupervised')),'vectorized_unsupervised' ); 
      xlswrite(fullfile(pwd,'\save\vectorized_mat_UNsupervised.xlsx'), vectorized_unsupervised, 'Sheet1', 'A1'); 
    end 
    
    function vectorize_unsupervised_pawel()
        load(fullfile(pwd, CreateTree.save_directory('swc_vector.mat')));
        load(fullfile(pwd, CreateTree.save_directory('vect_branching_pattern.mat')));
        load(fullfile(pwd, CreateTree.save_directory('vect_tortuosity.mat')));                
        load(fullfile(pwd, CreateTree.save_directory('vect_flux.mat')));               
        load(fullfile(pwd, CreateTree.save_directory('vect_leaf_index.mat')));
        load(fullfile(pwd, CreateTree.save_directory('vect_energy.mat')));
        load(fullfile(pwd, CreateTree.save_directory('vect_wiring.mat'))); 
        load(fullfile(pwd, CreateTree.save_directory('vect_taper_rate.mat')));

        vect_mat = {vect_energy vect_flux  vect_leaf_index  vect_branching_pattern  vect_wiring  vect_tortuosity vect_taper_rate  };
        parts = 2;
        vectorized_mat = zeros(length(swc_vector), parts*length(vect_mat));
        for j=1:length (vect_mat)
          grand_mat =[];
          mat = vect_mat(j);
           for k=1:length(mat{1})
               mat1 = cell2mat(mat{1}(k));
               grand_mat= [grand_mat ; mat1]
           end 
           
      
           min_val = min(grand_mat(:,1));
           max_val = max(grand_mat(:,1)); 
           
           % get the interval           
           val = (max_val- min_val)/parts;           
           interval = zeros(parts,1);
               for p =1:parts
                  interval(p,1) = p*val;
               end
       %loop through the descriptors        
     for e=1:length (vect_mat)
       mat = vect_mat(e);
          % compute the value of D for each T at interval points
          %loop through the trees in each descriptor
          for k=1:length(mat{1})
              col = (e*parts)-1;
               % loop through the trees
               mat1 = cell2mat(mat{1}(k));
         % for each tree get the values of teh function at the intervals
               mat1 = sortrows(mat1,1);
               %loop through the inetrvals 
               for a=1:length(interval) 
                   found =0;
                   for s=1:length(mat1)-1
                      if mat1(s,1)== interval(a) || mat1(s+1,1)== interval(a) 
                        vectorized_mat(k,col) = mat1(s,2);
                        col=col+1;
                        found =1;
                        break;
                      elseif mat1(s,1)<interval(a) && mat1(s+1,1)>interval(a)
                         vectorized_mat(k,col) =  mat1(s+1,2);
                         col=col+1;
                         found =1;
                         break; 
                      end 
                   end 
                   % if the interval value is larger that the largest
                   % radius then take the largest radius
                   if found ==0
                   vectorized_mat(k,col) =  mat1(end,2); 
                   col=col+1;
                   end 
               end 
          end
     end 
 end        
      save(fullfile( pwd , CreateTree.save_directory('vectorized_mat_pawel')),'vectorized_mat' ); 
      xlswrite(fullfile(pwd,'\save\vectorized_mat_UNsupervised_pawel.xlsx'), vectorized_mat, 'Sheet1', 'A1'); 
end 

    function combine_matrices_with_weights(weights,dist_mat)
%         %Now we create weights with which they should be used. 
%         weights = [1 1 1 1 1];
        
     combined_dist_mat = zeros( size(dist_mat(1)) );
        for i = 1:length( weights )
            mat = cell2mat(dist_mat(i));
            combined_dist_mat = combined_dist_mat + weights(i).*mat;
        end
        get_dendrogram(combined_dist_mat,'combined_dist_mat');
%         DM.get_clusters(combined_dist_mat)

    end 
    
    function new_counter = increment_counter( current , start, stop, dx )
        %Find a position where we can increment
        position_to_increment = 1;
        new_counter = current;
        while ( (position_to_increment ~= length(start)+1) && (new_counter(position_to_increment) >= stop(position_to_increment)) )
            position_to_increment = position_to_increment + 1;
        end

        %And, if there is still something to increment:
        if ( position_to_increment ~= length(start)+1 )        
            new_counter( position_to_increment ) = new_counter( position_to_increment )+dx;
           for kk = 1:(position_to_increment-1)
               new_counter(kk) = start(kk);                
           end
        end
    end
    
    
        
    end
end

