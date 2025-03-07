function B_values = BBs(S)
global m_exp;
global u_star;
global phi_ustar;

%phi_values = @(S) (S.^m_exp)./ ((1-S).^m_exp);

    Nn = length(S);
    B_values = zeros(1, Nn);

    for i = 1:Nn
         if S(i) < 0
            B_values(i) = S(i);
        elseif S(i) < u_star
            B_values(i) = max((S(i).^m_exp) ./ ((1 - S(i)).^m_exp), 0);
        else
            B_values(i) = phi_ustar + S(i) - u_star;
        end
    end
end