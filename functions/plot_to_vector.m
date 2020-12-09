function vect = plot_to_vector( x , y , min_x , max_x , numeber_of_steps )
%%%In this function we take a partially constant-plot, given by (x,y) coordinates along
%%%with the range in which the vectorization is happening and the number of
%%%steps of vectorization (the lenght of the output vector). This procedure
%%%simply iterate throught the range in between min_x and max_x and get the
%%%values of the given partially constant function therein. 
 
    %First we initialize the vector:
    vect = zeros(1,numeber_of_steps+1);
    for i = 1:(numeber_of_steps+1)
        vect(i) = y( length(y) );
    end
    
    dx = (max_x-min_x)/numeber_of_steps;
    xx = min_x;
    current_x = 1;    
    for i = 1:(numeber_of_steps+1)
        while ( (current_x < length(x)) && (x(current_x) < xx) )
           current_x = current_x + 1; 
        end
        %disp( current_x )
        if ( current_x == length(x) )
            %in this case, the value of the function will not change, so we
            %just leave the initial values put into the vect:
            break;
        end
        %in this case the value of the function will be y( current_x )
        %disp( y( current_x ) );
        vect(i) = y( current_x );       
        
        xx = xx + dx;
    end
end

%Test:
%x = [1,2,3]
%y = [3 2 1]
%return 3 3 2 1

%x = [0.5 1.5 2.5 3.5 4.5 , 5.5]
%y = [0 3 0.5 2 1 -1]
%return ???
