function B = PBBs(S, m)
% PBBs: Calculate B values based on S and parameter m with specific conditions.
% S: a 2D array of size nx+1 x ny+1
% m: a scalar value
% Returns:
% B: a 2D array of size nx+1 x ny+1 containing the computed B values.

[n, dim] = size(S); % Get the size of the input matrix

B = zeros(n, dim); % Initialize the output matrix
T = (1/m)^(1/(m-1));

for i = 1:n
    for j = 1:dim
        if S(i,j) < 0
            B(i,j) = S(i,j); % If S < 0, B = S
        elseif S(i,j) <= T
            B(i,j) = (S(i,j)).^m; 
        else
            B(i,j) = S(i,j) - T + T^m; 
        end
    end
end
end
