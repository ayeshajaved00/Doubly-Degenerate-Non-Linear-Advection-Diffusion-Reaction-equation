function [S] = InitialWW(X,m,t,d)
R=sqrt(X.^2);
alpha=1/(m-1 + 2/d);
gamma=1.5;

S=(1+t)^(-alpha)*(max(gamma - alpha*(m-1)*R.^2/(2*d*m*(t+1)^(2*alpha/d)),0)).^(1/(m-1));
end 