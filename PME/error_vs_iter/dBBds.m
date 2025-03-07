function Bp = dBBds(S,m)
n=length(S);
Bp=zeros(1,n);
T = (1/m)^(1/(m-1)); 

for i = 1:n
        if S(i) < 0
            Bp(i) = 1;
        elseif S(i) < T
            Bp(i) = m * S(i).^(m-1); 
        else
            Bp(i) = 1; 
        end
end
end
