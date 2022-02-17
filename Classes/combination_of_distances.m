classdef combination_of_distances  
   
methods(Static)
    function combine_distances(dx)
         tStart = tic; 
        %In this file we will read distance matrices among various descriptor
        %and build a combinations of them to get one distance matrix to rule them
        %all.
        %load swc vector
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
      
        tortuosity      = tortuosity_dist_matrix;        
        branching       = branching_dist_matrix;        
        flux            = flux_dist_matrix;      
        leaf_index      = leaf_index_dist_matrix;        
        energy          = energy_dist_matrix;
        wiring          = wiring_dist_matrix;
        taper      = taper_rate_dist_matrix;
         tmd      = tmd_dist_matrix;
        
       %Now we create weights with which they should be used. 
        dist_mat_1 = {wiring energy tortuosity branching flux leaf_index tmd taper};
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
        % %Below a prototype of code to get optimal weights.        
        start = zeros( 1,length(dist_mat) );
        current = start;
        stop = repmat(3,1,length(dist_mat_1))

        current_best_score = 0;
        current_best_weights = start;
		
        %swc_vector is packed with info             
        %classes: has the unique values . i.e. the number of neuron classes
        %idx: is an array of indices of the occurances in these classes
		[classes,~,idx]=unique(swc_vector(4,:)');
        while ( 1 )
            %We will use the counter current as weights: 
            %Construc the combined distance matrix:
            combined_dist_mat = zeros( size(dist_mat(1)) );
                for i = 1:length(current)
                    mat = cell2mat(dist_mat(i));
                    combined_dist_mat = combined_dist_mat + current(i).*mat;
                end
                   
            %get maximal distance within each class
            mat=[];
                for i=1:length(classes)
                    % get all the indices of the class i 
                    k = find(idx==i);
                    %starting index
                    s = min(k);
                    %ending index
                    e = max(k); 
                    mat = [ mat; [i s e]];
                    max_distance_within_i_class(i) = max( max( combined_dist_mat(s:e,s:e)));
                end 
            internal_dist =  sum(max_distance_within_i_class);
            
            %get the external distance
            % choose pairs 
            c = nchoosek(mat(:,1),2);
            [r,~]=size(c);
                for i=1:r
                    % start and end of pair 1
                    s1=mat(c(i,1),2);
                    e1=mat(c(i,1),3);
                    % start and end of pair 2
                    s2=mat(c(i,2),2);
                    e2=mat(c(i,2),3);
                    min_distance_between_i_j_classes(i) =min(min(combined_dist_mat(s1:e1,s2:e2)));        
                end    
            external_dist = min(min_distance_between_i_j_classes);

            score = external_dist / internal_dist;

            if ( score >= current_best_score )
                current_best_score = score;
                current_best_weights = current;
            end

            new_counter = combination_of_distances.increment_counter( current , start, stop, dx );
            if ( current == new_counter )
                break;
            end
            current = new_counter;
            disp(current_best_weights);
            disp(score);
        end 
        assignin('base','score',score);
        assignin('base','current_best_weights',current_best_weights);
        assignin('base','combined_dist_mat',combined_dist_mat);
        save(fullfile( pwd , CreateTree.save_directory('score')),'score' );
        save(fullfile( pwd , CreateTree.save_directory('current_best_weights')),'current_best_weights' );
        save(fullfile( pwd , CreateTree.save_directory('combined_dist_mat')),'combined_dist_mat' );
        combination_of_distances.combine_matrices_with_weights(current_best_weights,dist_mat); 
        
        tEnd = toc(tStart);
        fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
    end


    function combine_matrices_with_weights(weights,dist_mat)
        
     combined_dist_mat = zeros( size(dist_mat(1)) );
        for i = 1:length( weights )
            mat = cell2mat(dist_mat(i));
            combined_dist_mat = combined_dist_mat + weights(i).*mat;
        end
        get_dendrogram(combined_dist_mat,'combined_dist_mat');

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

