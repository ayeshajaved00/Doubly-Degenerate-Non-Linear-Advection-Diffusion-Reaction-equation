function b = Pbs(S, m)

[n, dim] = size(S); % Get the size of the input matrix

b = zeros(n, dim); % Initialize the output matrix
T = (1/m)^(1/(m-1)); % Threshold value based on m

for i = 1:n
    for j = 1:dim
        if S(i,j) <= T
            b(i,j) = S(i,j); 
        else
            b(i,j) = (S(i,j) - T + T^m).^(1/m); 
        end
    end
end

end
