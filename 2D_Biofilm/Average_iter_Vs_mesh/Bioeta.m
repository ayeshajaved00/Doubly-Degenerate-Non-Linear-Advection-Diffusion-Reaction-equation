function [Y]=Bioeta(E,a,b,dA)

Q=sqrt(sum(sum((E.^(1/2).*a+E.^(-1/2).*b).^2*dA)));

R=sqrt(sum(sum((-E.^(1/2).*a+E.^(-1/2).*b).^2*dA)));

Y = (Q+R)/2;
end
