function BB_values = BiodBBds(S)
    global m_exp;
    global u_star;
    global phi_ustar;

    [n, dim] = size(S);
    BB_values = zeros(n, dim);

    for i = 1:n
        for j = 1:dim
            s = S(i, j);
            if s < 0
                BB_values(i, j) = 1;
            elseif s < u_star
                BB_values(i, j) = max(m_exp * (s^(m_exp - 1)) / (1 - s)^(m_exp + 1), 0);
            else
                BB_values(i, j) = 1;
            end
        end
    end
end
