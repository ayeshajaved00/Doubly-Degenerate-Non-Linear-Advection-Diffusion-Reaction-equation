function bb_values = Biodbds(S)
    global m_exp;
    global u_star;
    global phi_ustar;

    [n, dim] = size(S);
    bb_values = zeros(n, dim);

    for i = 1:n
        for j = 1:dim
            if S(i,j) < 0
                bb_values(i,j) = 0; 
            elseif S(i,j) < u_star
                bb_values(i,j) = 1;
            else
                tmp = (phi_ustar + S(i,j) - u_star).^(1/m_exp);
                bb_values(i,j) = 1./(m_exp * tmp.^(m_exp - 1) .* (1 + tmp).^2) + 1e-10; % Avoid division by zero
            end
        end
    end
end
