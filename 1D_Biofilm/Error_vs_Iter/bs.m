function b_values = bs(S)
global m_exp;
global u_star;
global phi_ustar;

    Nn = length(S);
    b_values = zeros(1,Nn);

    for i = 1:Nn
         if S(i) < 0
            b_values(i) = 0;
        elseif S(i) < u_star
            b_values(i) = max(S(i), 0); 
        else
            tmp = max(((phi_ustar + S(i) - u_star)).^(1/m_exp), 0);
            b_values(i) = tmp ./ (1 + tmp);
        end
    end
end