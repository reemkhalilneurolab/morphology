function dist = distance_of_vectors( v1 , v2 )
    if ( length(v1) ~= length(v2) )
        msg = 'Error, vectors of different sizes in compute distance procedire.';
        error(msg)
    end
    dist = 0;
    for i = 1:length(v1)
        dist = dist + (v1(i)-v2(i))^2;
    end
    dist = sqrt(dist);
end