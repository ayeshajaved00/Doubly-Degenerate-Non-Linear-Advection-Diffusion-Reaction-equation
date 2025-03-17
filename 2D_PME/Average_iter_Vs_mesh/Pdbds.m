function bpp = Pdbds(S, m)

[n, dim] = size(S); % Get the size of the input matrix

bpp = zeros(n, dim); % Initialize the output matrix
T = (1/m)^(1/(m-1)); 

for i = 1:n
    for j = 1:dim
        if S(i,j) <= T
            bpp(i,j) = 1; 
        else
            bpp(i,j) = (1/m) .* (S(i,j) - T + T^m).^(1/m - 1); 
        end
    end
end

end
