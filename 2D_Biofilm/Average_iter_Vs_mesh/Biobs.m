function b_values = Biobs(S)
global m_exp;
global u_star;
global phi_ustar;

[n, dim] = size(S);

b_values = zeros(n, dim);

for i = 1:n
    for j = 1:dim
        if S(i,j) < 0
            b_values(i,j) = 0;
        elseif S(i,j) < u_star
            b_values(i,j) = max(S(i,j), 0); 
        else
            tmp = max(((phi_ustar + S(i,j) - u_star)).^(1/m_exp), 0);
            b_values(i,j) = tmp ./ (1 + tmp);
        end
    end
end
end
