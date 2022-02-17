function [ m ,b ] = get_linear_regression(x,y)
    P=polyfit(x,y,1);
    m=P(1);
    b=P(2);
    % linefit = polyval(P,x)
end 