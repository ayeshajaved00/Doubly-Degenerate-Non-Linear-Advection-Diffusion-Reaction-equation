function Fs_derivative = dFds(b_iter, dbds_iter)
    global lambda

    % Ensure b_iter values are not too small to avoid numerical issues
    b_iter = max(b_iter, 1e-9);

    % Compute intermediate terms
    T = (1 - b_iter.^(1 / lambda)).^lambda;
    Y = (1 - b_iter.^(1 / lambda)).^(lambda - 1);

    % Initialize derivative array
    Fs_derivative = zeros(1, length(b_iter));

    % Compute Fs_derivative for each element
    for i = 1:length(b_iter)
        term1 = 2 * (1 - T(i)) * dbds_iter(i) * (b_iter(i))^((-lambda + 1) / lambda) * sqrt(b_iter(i)) * Y(i);
        term2 = (dbds_iter(i) * (1 - T(i))^2) / (2 * sqrt(b_iter(i)));
        
        Fs_derivative(i) = term1 + term2;
    end
end
