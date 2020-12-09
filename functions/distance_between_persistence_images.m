 function dist = distance_between_persistence_images( image1 , image2 , p )
    %--------------------------------------
    %add comment
    %-----------------
    a = size(image1);
    b = size(image2);
    if ( a ~= b )
        throw "Non compatible sizes of images.";
    end
    dist = 0;
    for i=1:a(1)
        for j=1:a(2)
            dist = dist + ( image1(i,j)-image2(i,j) )^p;
        end
    end
end