function b = Pbs(S, m)
% Pbs: Compute the function B(S) based on the input S and parameter m.
% S: a 2D array of size nx+1 x ny+1
% m: a scalar value
% Returns:
% b: a 2D array of the same size as S with the computed values.

[n, dim] = size(S); % Get the size of the input matrix

b = zeros(n, dim); % Initialize the output matrix
T = (1/m)^(1/(m-1)); % Threshold value based on m

for i = 1:n
    for j = 1:dim
        if S(i,j) < 0
            b(i,j) = 0; % If S < 0, B(S) = 0
        elseif S(i,j) <= T
            b(i,j) = S(i,j); 
        else
            b(i,j) = (S(i,j) - T + T^m).^(1/m); 
        end
    end
end

end
