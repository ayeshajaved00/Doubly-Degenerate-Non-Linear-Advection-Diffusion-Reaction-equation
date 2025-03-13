function Bp = dBBds(S)
    n = length(S);
    Bp = zeros(1, n);

    for i = 1:n
        if S(i) <= 0
            Bp(i) = 0;  % Derivative is zero for S <= 0
        elseif (S(i) > 0) && (S(i) <= 1/sqrt(2))
            Bp(i) = S(i) ./ sqrt(1 - S(i).^2);  % Derivative in region 0 < S <= 1/sqrt(2)
        else
            Bp(i) = 1;  % Derivative is 1 for S > 1/sqrt(2)
        end
    end
end

