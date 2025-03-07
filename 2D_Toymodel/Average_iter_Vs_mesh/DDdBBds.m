function BBp = DDdBBds(S)

[n, dim] = size(S); 

BBp = zeros(n, dim);
for i = 1:n
    for j = 1:dim
        if S(i,j) < 0
            BBp(i,j) = 1;
        elseif (S(i,j) > 0) && (S(i,j) < 1/sqrt(2))
            BBp(i,j) = S(i,j)./(sqrt(1-S(i,j).^2));
        else
            BBp(i,j) = 1;
        end
    end
end
end

