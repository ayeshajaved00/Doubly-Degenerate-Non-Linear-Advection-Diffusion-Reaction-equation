function [Y]=eta(E,a,b,dx)

Q=sqrt(sum(sum(((E.^(1/2).*a+E.^(-1/2).*b).^2)*dx)));

R=sqrt(sum(sum(((-E.^(1/2).*a+E.^(-1/2).*b).^2)*dx)));

 Y = (Q+R)/2;

end
