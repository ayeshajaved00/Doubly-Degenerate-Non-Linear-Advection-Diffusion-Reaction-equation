function B_values = BioBBs(S)
global m_exp;
global u_star;
global phi_ustar;

[n, dim] = size(S);

B_values = zeros(n, dim);

for i = 1:n
    for j = 1:dim
        if S(i,j) < 0
            B_values(i,j) = S(i,j);
        elseif S(i,j) < u_star
            B_values(i,j) = max(((S(i,j)).^m_exp) ./ ((1 - S(i,j)).^m_exp), 0);
        else
            B_values(i,j) = phi_ustar + S(i,j) - u_star;
        end
    end
end
end
