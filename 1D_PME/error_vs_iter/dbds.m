function b = dbds(S,m)

n=length(S);
b=zeros(1,n);

T = (1/m)^(1/(m-1)); 

for i = 1:n
        if S(i) <= T
            b(i) = 1; 
        else
            b(i) = (1/m) .* (S(i) - T + T^m).^(1/m - 1); 
        end
end
end

  
