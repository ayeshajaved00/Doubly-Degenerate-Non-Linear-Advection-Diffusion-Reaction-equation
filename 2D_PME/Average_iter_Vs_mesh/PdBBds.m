function BBp = PdBBds(S, m)

[n, dim] = size(S); % Get the size of the input matrix

BBp = zeros(n, dim); % Initialize the output matrix
T = (1/m)^(1/(m-1)); % Threshold value based on m

for i = 1:n
    for j = 1:dim
        if S(i,j) <= T
            BBp(i,j) = m .* S(i,j).^(m-1); 
        else
            BBp(i,j) = 1; 
        end
    end
end
end
