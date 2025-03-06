clear all
close all
 clc
 
 %% Parameter and function definitions

%van Genuchten parameter
 
global Ks
global Sp
global lambda

lambda = 0.8;

%Permeability and saturation function

Ks = @(u, lambda) (sqrt(u).*(1-(1-u.^(1/lambda)).^lambda).^2);
Bp = @(u,p, lambda) min(Ks(u,lambda)/Sp(p, lambda),1);
bp = @(u,p, lambda) min(Sp(p, lambda)/Ks(u,lambda),1);

% Defining global arrays of u-values

global du
global u_star
global u_start

du = 1e-5;
u_start=1e-5;
u_values =u_start:du:1;

% Calculate pressure for each u
p_values = 1 - ((u_values.^(-1/lambda) - 1).^(1 - lambda));

%% Find corresponding $w=\Phi(u)$-values 

w_values(1)=0;
b_values(1)=0;
B_values(1)=0;

%Integrate over the subintervals [0, p(i)]

for i = 2:length(p_values)
    % For the subsequent subintervals [p(i-1), p(i)]
    hp = (p_values(i) - p_values(i-1));
    w_values(i) = w_values(i-1) + (hp/8) *((Ks(u_values(i-1), lambda) + 3*Ks((2*u_values(i-1)+u_values(i))/3 , lambda) + 3*Ks((u_values(i-1)+2*u_values(i))/3,lambda) + Ks(u_values(i), lambda)));   
end

%Find $\Phi'(u)$ from the w-values
Phip = gradient(w_values,u_values);

%Find $u^*$ and $\Phi(u^*)$
u_star = min(find(Phip>1));
phi_ustar = w_values(u_star);
u_star = u_values(u_star);

%% calculate $b(s)$ and $B(s)$ corresponding to $s$ values
global s
global s_end
global N
%The s-input string 

s_end=2;
s = u_start:du:s_end;
N = length(s);
b_values = zeros(1,N);
B_values = zeros(1,N);

%Computing $B$-values corresponding to the s-values by using $B(s)=s$ for
%$s<u^*$ and $B(s)=\Phi(u^*)+ (s-u^*)$
for i = 1:N
    if s(i) < u_star
        B_values(i) = w_values(i);
    else
        B_values(i) = phi_ustar + s(i) - u_star;
    end
end

% Computing the  $b$-values by linearly interpolating between the w-values
j=2;
w_at_1=w_values(end);

for i = 1:N
    
    while (w_values(j)<B_values(i))&&(B_values(i)<w_at_1)
             j=j+1;
    end

    b_values(i)= (u_values(j) * ( B_values(i)- w_values(j-1)) + u_values(j-1) * (w_values(j)-B_values(i))) / (w_values(j) - w_values(j-1));
end
