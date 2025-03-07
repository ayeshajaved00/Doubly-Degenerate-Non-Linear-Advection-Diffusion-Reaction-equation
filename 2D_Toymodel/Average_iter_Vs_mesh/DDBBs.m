function B = DDBBs(S)
% S: a 2D array of size nx+1 x ny+1
% m: a scalar value
% nx: number of elements in x-direction
% ny: number of elements in y-direction
% b: a 2D array of size nx+1 x ny+1, containing the values calculated by the given code function
%B=max(0,S.^m);
[n, dim] = size(S); % Get the size of the input vector

B = zeros(n, dim);

for i = 1:n
    for j = 1:dim
        if S(i,j) <= 0
            B(i,j) = S(i,j); %before B(i,j) =0
        elseif (S(i,j) >= 0) && (S(i,j) <= 1/sqrt(2))
            B(i,j) =  1-sqrt(1-S(i,j).^2);
         else
             B(i,j) = S(i,j)+1-sqrt(2);
        end
    end
end
