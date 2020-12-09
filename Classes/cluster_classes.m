classdef cluster_classes
   
methods(Static)
    function cluster
        tStart = tic; 
        %In this file we will read distance matrices comong from various descriptor
        %and build a combinations of them to get one distance matrix to rule them
        %all.
        load_switch =1;
%         load swc vector
        if load_switch ==1
            load(fullfile(pwd, CreateTree.save_directory('swc_vector.mat')));
            %First we read the distance matrices:
            load(fullfile(pwd, CreateTree.save_directory('branching_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('tortuosity_dist_matrix.mat')));                
            load(fullfile(pwd, CreateTree.save_directory('flux_dist_matrix.mat')));               
            load(fullfile(pwd, CreateTree.save_directory('leaf_index_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('energy_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('wiring_dist_matrix.mat')));
            load(fullfile(pwd, CreateTree.save_directory('taper_rate_dist_matrix.mat')));
             load(fullfile(pwd, CreateTree.save_directory('tmd_dist_matrix.mat')));
    end 
      
        tortuosity      = tortuosity_dist_matrix;        
        branching       = branching_dist_matrix;        
        flux            = flux_dist_matrix;      
        leaf_index      = leaf_index_dist_matrix;        
        energy          = energy_dist_matrix;
        wiring          = wiring_dist_matrix;
        taper      = taper_rate_dist_matrix;
         tmd      = tmd_dist_matrix;  

         dist_mat = {energy flux leaf_index branching wiring tortuosity tmd taper};
        [classes,~,idx]=unique(swc_vector(4,:)');
        mat=[];
                for i=1:length(classes) %neurons classes
                    % get all the indices of the class i 
                    k = find(idx==i);
                    %starting index
                    s = min(k);
                    %ending index
                    e = max(k); 
                    mat = [ mat; [i s e]];
                end              

        [r1 c1]=size(dist_mat{1})

        for r=1:r1
          for c=1:c1
            vect =[dist_mat{1}(r,c) dist_mat{2}(r,c) dist_mat{3}(r,c) dist_mat{4}(r,c) dist_mat{5}(r,c) ...
                   dist_mat{6}(r,c) dist_mat{7}(r,c) dist_mat{8}(r,c)];
                   
               
               mini(r,c)= min(vect);              
          end 
            
        end 

         dist_mat{end+1}= mini;

                [rw col]=size(mat);
         final_detection_rate =[];
         for j=1:length(dist_mat) % loop through the descriptors 
%            
            for v=1:rw  % loop through the classes within each descriptor
                %number of Neurons in the class
                rate_of_detection_i=[];
                for b=mat(v,2):mat(v,3)
                    proportions=[];
                    internal_vect= dist_mat{j}(b,mat(v,2):mat(v,3));
                    external_vect=[];
                       for t=1:rw
                           if t~=v
                           external_vect= [external_vect dist_mat{j}(b,mat(t,2):mat(t,3))];
                           end 
                       end 
                   internal_vect = sort(internal_vect,'descend');
                   [r,c] = size(internal_vect);
                   [r1,c1] = size(external_vect);
                   for k=1:c
                        eps= internal_vect(k); 
                            internal_count = sum(internal_vect<=eps);
                            internal_ratio = internal_count/c; % m/m = 100% (all neurons in C_1 within distance d_m from N_1)
                            external_count = sum(external_vect<=eps);
                            external_ratio = internal_count/(internal_count+external_count); % m/c_m+m (the percentage of all neurons that are a distance d_m from N_1).
                            proportions = [proportions min(internal_ratio,external_ratio)];  % the smaller of the two is  beta_m = m/c_m+m
                            %At the end of the iterative process (associated to N_1) we get the proportions beta_1, beta_2,.... , beta_m
%                           end 
                   end 
                  rate_of_detection_i(b)= max(proportions);
                end  %end  for b=mat(v,3):mat(v,2)
                rate_of_detection_all(v)= max(rate_of_detection_i);
            end % for v=1:rw
          final_detection_rate= [final_detection_rate ; rate_of_detection_all];
         end %j=1:length(dist_mat) 
          
       assignin('base','final_detection_rate',final_detection_rate);   
        save(fullfile( pwd , CreateTree.save_directory('final_detection_rate')),'final_detection_rate' );
       xlswrite(fullfile(pwd,'\save\final_detection_rate.xlsx'), final_detection_rate, 'Sheet1', 'A1');  
            
       
        
    end 
    
end

end
