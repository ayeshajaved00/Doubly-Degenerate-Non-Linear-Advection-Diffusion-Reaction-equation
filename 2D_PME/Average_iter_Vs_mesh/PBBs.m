function B = PBBs(S, m)

[n, dim] = size(S); % Get the size of the input matrix

B = zeros(n, dim); % Initialize the output matrix
T = (1/m)^(1/(m-1));

for i = 1:n
    for j = 1:dim
        if S(i,j) <= T
            B(i,j) = (S(i,j)).^m; 
        else
            B(i,j) = S(i,j) - T + T^m; 
        end
    end
end
end
