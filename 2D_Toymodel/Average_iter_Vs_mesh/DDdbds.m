function bpp = DDdbds(S)

[n, dim] = size(S); 

bpp = zeros(n, dim);

% Loop through each element in S to calculate bpp
for i = 1:n
    for j = 1:dim
        if S(i, j) <= 0
            bpp(i, j) = 0; 
        elseif S(i, j) < 1/sqrt(2)
            bpp(i, j) = 1;
        elseif (S(i, j) > 1/sqrt(2)) && (S(i, j) < sqrt(2))
            bpp(i, j) = (sqrt(2) - S(i, j)) / sqrt(1 - (sqrt(2) - S(i, j))^2);
        else
            bpp(i, j) = 0;
        end
    end
end

end
