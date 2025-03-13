function B = BBs(S)
    n = length(S);
    B = zeros(1, n);

    for i = 1:n
        if S(i) <= 0
            B(i) = 0; 
        elseif (S(i) >= 0) && (S(i) <= 1/sqrt(2))
            B(i) = 1 - sqrt(1 - S(i).^2);
        else
            B(i) = S(i) + 1 - sqrt(2);
        end
    end
end


   
