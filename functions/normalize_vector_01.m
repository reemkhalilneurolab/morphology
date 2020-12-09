function normalized_vector = normalize_vector_01( vect )
    normalized_vector = zeros( 1,length(vect) );
    min_= min(vect);
    max_= max(vect);
    for i=1:length( vect )
        normalized_vector(i) = (vect(i)-min_)/( max_ - min_ );
    end
end