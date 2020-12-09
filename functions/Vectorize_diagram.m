 function result = Vectorize_diagram( diagram , min_x , max_x , min_y , max_y , size_x , size_y )
    %--------------------------------------
    %add comment
    %-----------------
    a = size(diagram);
    result = zeros( size_x , size_y );
    dx = ( max_x - min_x )/size_x;
    dy = ( max_y - min_y )/size_y;
    for xit = 1:size_x
        x = min_x + xit*dx;
        for yit = 1:size_y
            y = min_y + yit*dy;
            pt = [x y];
            %compute value at the point (x,y):
            sum = 0;
            for i = 1:a(1)
                %compute distance between (x,y) and ( diagram(i,1),diagram(i,2) )
                d = distance_of_vectors(pt ,diagram(i,:) );
                sum = sum + 1/(d^2+1);                
            end
            if sum == 0
                disp("zero");
            end
            result( xit,yit ) = sum;
        end
    end
end