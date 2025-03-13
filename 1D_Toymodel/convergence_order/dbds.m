function b = dbds(S)

n=length(S);
b=zeros(1,n);

for i = 1:n
        if S(i) < 1/sqrt(2)
            b(i) = 1;
        elseif (S(i) >= 1/sqrt(2)) && (S(i) <= sqrt(2))
            z = sqrt(2) - S(i); b(i) = z ./ sqrt(1 - z.^2);
        else
            b(i) = 0;
        end
end
end
