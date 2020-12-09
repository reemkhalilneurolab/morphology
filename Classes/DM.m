classdef DM
    
    methods(Static)
%%
%%        
        function distance = get_distance_of_plots_grid_version( x1,y1,x2,y2 )
          A_class = zeros(length(x1),1);
          B_class = ones(length(x2),1);
          A = horzcat(x1,y1,A_class);
          B = horzcat(x2,y2,B_class);
          C = vertcat(A,B);
          C= sortrows(C,1);
          distance =0;
          
          for i=1:size(C,1)-1
           
              d= C(i+1,1)-C(i,1);
              xi_class =  C(i,3);
              yi= C(i,2);
              cnt =i;
              if cnt >1 
                 while (cnt >1)  % loop down for all the values lower than xi
                   cnt=cnt-1;
                   if C(cnt,3) ~= xi_class
                     y_first= C(cnt,2); % get teh y  
                     break;
                   else
                    y_first =0;
                   end 
                 end  
              else               
                  y_first =0;
              end               
            distance = distance + d*abs(yi-y_first) ;
          end 
          
        end 
%%
%%       
        function distance = get_add_distance_of_plots_grid_version( x1,y1,x2,y2 )
        distance = 0;
        previous_x = 0;
        x1_counter = 1;
        x2_counter = 1;
        while ( (x1_counter<=length(x1)) && (x2_counter<=length(x2)) )
            if ( x1(x1_counter) < x2(x2_counter) ) 
                %if x1<x2 subtract the y's abs(y2-y1) and multiply by x1 (
                %u get the area of the rectangle under the step fuction
                %beween x1 and x2
                distance = distance + abs(y1(x1_counter)+y2(x2_counter))*(x1(x1_counter)-previous_x);
                previous_x = x1(x1_counter);
                x1_counter = x1_counter+1;
            else           
                distance = distance + abs(y1(x1_counter)+y2(x2_counter))* (x2(x2_counter)-previous_x);
                previous_x = x2(x2_counter);
                x2_counter = x2_counter+1;
            end        
        end
        while ( x1_counter<=length(x1) )
            distance = distance + abs(y1(x1_counter))*(x1(x1_counter)-previous_x);
            previous_x = x1(x1_counter);
            x1_counter = x1_counter+1;
        end
        while ( x2_counter<=length(x2) )
            distance = distance + abs(y2(x2_counter))*(x2(x2_counter)-previous_x);
            previous_x = x2(x2_counter);
            x2_counter = x2_counter+1;        
        end

    end%get_distance_of_plots_grid_version
%%
%%    
    function get_heat_map(distance_matrix)         
        [r c] = size(distance_matrix);
        H1=figure;
        imagesc(distance_matrix);
        ylabel('set1');
        xlabel('set2');
        HC = colorbar;
        ylabel(HC,'distance');         
    end 
%%
 
    end
end



