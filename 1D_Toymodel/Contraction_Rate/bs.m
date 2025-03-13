function b = bs(S)
    n = length(S);
    b = zeros(1, n);

    for i = 1:n
        if S(i) <= 1/sqrt(2)
            b(i) = S(i);
        elseif (S(i) >= 1/sqrt(2)) && (S(i) <= sqrt(2))
            b(i) = sqrt(1 - (sqrt(2) - S(i)).^2);
        else
            b(i) = 1;
        end
    end
end
