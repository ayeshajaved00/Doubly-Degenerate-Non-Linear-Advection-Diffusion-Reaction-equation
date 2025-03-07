function bpp = Pdbds(S, m)
% Pdbds: Compute the derivative of B(S) with respect to S based on the input S and parameter m.
% S: a 2D array of size nx+1 x ny+1
% m: a scalar value
% Returns:
% bpp: a 2D array of the same size as S with the computed derivative values.

[n, dim] = size(S); % Get the size of the input matrix

bpp = zeros(n, dim); % Initialize the output matrix
T = (1/m)^(1/(m-1)); 

for i = 1:n
    for j = 1:dim
        if S(i,j) < 0
            bpp(i,j) = 0; % If S < 0, b'(S) = 0
        elseif S(i,j) < T
            bpp(i,j) = 1; 
        else
            bpp(i,j) = (1/m) * (S(i,j) - T + T^m)^(1/m - 1); 
        end
    end
end

end
