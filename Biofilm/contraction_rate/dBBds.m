 function BB_values = dBBds(S)
global m_exp;
global u_star;
global phi_ustar;

%phi_prime_values = @(S) (m_exp*S.^(m_exp-1))./ ((1-S).^(m_exp+1));

%u_star = 0.326220512390137;

    Nn = length(S);

    BB_values = zeros(1, Nn);

    for i = 1:Nn
            if S(i) < 0
                BB_values(i) = 1;
            elseif S(i) < u_star
                BB_values(i) = max((m_exp * (S(i).^(m_exp - 1))) ./ ((1 - S(i)).^(m_exp + 1)), 0);
            else
                BB_values(i) = 1;
            end
    end
end
 