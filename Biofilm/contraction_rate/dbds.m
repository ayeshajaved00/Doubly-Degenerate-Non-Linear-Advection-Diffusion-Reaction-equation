function bb_values = dbds(S)

global m_exp;
global u_star;
global phi_ustar;

%p = (phi_ustar + S - u_star).^(1/4);

    Nn = length(S);

    bb_values = zeros(1, Nn);

    for i = 1:Nn
           if S(i) < 0
                bb_values(i) = 0; % Condition for S(i) < 0
            elseif S(i) < u_star
                bb_values(i) = 1; % Condition for S(i) < u_star
            else
                tmp = (phi_ustar + S(i) - u_star).^(1/m_exp);
                bb_values(i) = 1./(m_exp .* tmp.^(m_exp - 1) .* (1 + tmp).^2)+ 1e-10; % Avoid division by zero
            end
    end
end
