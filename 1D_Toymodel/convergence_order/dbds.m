function b = dbds(S)

n=length(S);
b=zeros(1,n);

for i = 1:n
        if S(i) <= 0
            b(i) = 1e-10;
        elseif S(i) < 1/sqrt(2)
            b(i) = 1;
        elseif (S(i) > 1/sqrt(2)) && (S(i) < sqrt(2))
            b(i) = (sqrt(2) - S(i)) / sqrt(1 - (sqrt(2) - S(i))^2);
        else
            b(i) = 0;
        end
end
end
