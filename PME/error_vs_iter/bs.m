function b = bs(S,m)

n=length(S);
b=zeros(1,n);
T = (1/m)^(1/(m-1));

for i = 1:n
        if S(i) < 0
            b(i) = 0;
        elseif S(i) <= T
            b(i) = S(i); 
        else
            b(i) = (S(i) - T + T^m).^(1/m); 
        end
end
end



