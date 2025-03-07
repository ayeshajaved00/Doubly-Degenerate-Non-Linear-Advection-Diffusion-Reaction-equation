function [S] = DDBarenblatt(X,Y,m,t,d)
% Returns a 2D array of values of the InitialWW function at the given x and y
% values with the given parameters Plx, Prx, Ply, and m.
R=sqrt(X.^2+Y.^2);
alpha=1/(m-1 + 2/d);
gamma=1.5;

S=(1+t)^(-alpha)*(max(gamma - alpha*(m-1)*R.^2/(2*d*m*(t+1)^(2*alpha/d)),0)).^(1/(m-1));

end
