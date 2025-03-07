function BBp = PdBBds(S, m)
% PdBBds: Calculate the derivative of B(S) with respect to S.
% S: a 2D array of size nx+1 x ny+1
% m: a scalar value
% Returns:
% BBp: a 2D array of size nx+1 x ny+1 containing the computed derivative values.

[n, dim] = size(S); % Get the size of the input matrix

BBp = zeros(n, dim); % Initialize the output matrix
T = (1/m)^(1/(m-1)); % Threshold value based on m

for i = 1:n
    for j = 1:dim
        if S(i,j) < 0
            BBp(i,j) = 1; % If S < 0, B'(S) = 1
        elseif S(i,j) < T
            BBp(i,j) = m * S(i,j).^(m-1); 
        else
            BBp(i,j) = 1; 
        end
    end
end
end
