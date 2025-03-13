function b = DDbs(S)
[n, dim] = size(S); 

b = zeros(n, dim);

% Loop through each element in S to calculate b
for i = 1:n
    for j = 1:dim 
        if S(i, j) <= 1/sqrt(2)
            b(i, j) = S(i, j);
        elseif (S(i, j) >= 1/sqrt(2)) && (S(i, j) <= sqrt(2))
            b(i, j) = sqrt(1 - (sqrt(2) - S(i, j)).^2);
        else
            b(i, j) = 1;
        end
    end
end

end
