function the_same = compare_counters( first, second )
    the_same = 1;
    if ( length(first) ~= length(second) )
        the_same = 0;
    end
    for i = 1:length(first)
        if ( first(i) ~= second(i) )
            the_same = 0;
        end
    end
end