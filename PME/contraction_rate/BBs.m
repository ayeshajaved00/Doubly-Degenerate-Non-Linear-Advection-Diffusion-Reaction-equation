function B = BBs(S,m)
n=length(S);
B=zeros(1,n);
T = (1/m)^(1/(m-1));

for i = 1:n
        if S(i) < 0
            B(i) = S(i); 
        elseif S(i) <= T
            B(i) = (S(i)).^m; 
        else
            B(i) = S(i) - T + T^m; 
        end
end
end


   
