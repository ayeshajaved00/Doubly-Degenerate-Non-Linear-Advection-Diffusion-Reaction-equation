%% Biofilm Model - Contraction rate
% This script solves the Biofilm model using iterative 
% scheme (M scheme) and visualizes the contraction rate.

clc;
clear all;
close all;

%% Input values
global m_exp;
global u_star;
global phi_ustar;
m_exp = 6;
u_star = 0.367774963378906;
phi_ustar = 0.0387;

L = input('Left starting point of the domain: ');
R = input('Right starting point of the domain: ');
n = input('Number of cells: ');
T = input('Total runtime: ');
dt = input('Time step: ');

dx = (R - L) / n;

Mb = 0.01; 

MB = 0.01;

ee = 0.5;

x = zeros(1,n+1);
for i = 1:n+1
   x(i) = L+dx*(i-1); 
end

%% The Number of Time Steps

N = floor(T/dt)+1;  
t = zeros(1,N); 

%% Boundary condition
Pl = 0;
Pr = 0;

%% Initial condition
S0 = InitialWW(x, m_exp, 0, 1);
b_n = S0;
S_current = S0;

%% Basic iteration information
t_start = 1; t_end = N;
EA = [];
contraction_rates = []; 
errors = [];
iteration_no = 50;
diverge = 0;
%% Time loop
tot_iter_number = 0;
for j = t_start:t_end
    no_it_pertimestep = 0;
    t(j) = (j-1) * dt;

     if rem(j, 10) == 0
         contraction_rates = [];
         errors = [];
     end

    S_iter = S_current;

    %% Precomputing arrays
    tm = (dt / dx^2) * ones(1, n-1);
    tp = (dt / dx^2) * ones(1, n-1);

    Main_D = [1, tm + tp, 1];
    Lower_D = [0, -tm, 0];
    Upper_D = [0, -tp, 0];

   for i = 1:iteration_no
        b_iter = bs(S_iter);      % Value of b at S^{i-1}_n
        B_iter = BBs(S_iter);     % Value of B at S^{i-1}_n
        dbds_iter = dbds(S_iter); % Derivative b' at S^{i-1}_n
        dBds_iter = dBBds(S_iter);% Derivative B' at S^{i-1}_n

        %% Defining Lb and LB
        % M-scheme
        Lb = min(max(dbds_iter + Mb * dt, 2 * Mb * dt), 1);
        LB = min(max(dBds_iter + MB * dt, 2 * MB * dt), 1);

        % Matrix assembly
        AA = Main_D + (Lb ./ LB) .* (1 - ee * (1 - b_n) * dt);
        AA(1) = 1; 
        AA(end) = 1;

        f = b_n(2:n) - (1 - ee * (1 - b_n(2:n)) * dt) .* b_iter(2:n) ...
            + (1 - ee * (1 - b_n(2:n)) * dt) .* B_iter(2:n) .* (Lb(2:n) ./ LB(2:n));
        F = [Pl, f, Pr];
        if i == 1
           W_intiter = B_iter;
        else
           W_intiter = W_iter;
        end
        % Solve tridiagonal system
        W_iter = tridiagQQ(AA, Lower_D, Upper_D, F)';

        % Update S_iter
        S_intiter = S_iter;
        S_iter = (1 ./ LB) .* (W_iter - B_iter) + S_intiter;

        %% Error calculation
        dW_dx = gradient(W_iter - W_intiter, dx);

        err = sqrt(trapz(x, (S_iter - S_intiter).^2 + dt * (dW_dx.^2)));

        EA = [EA, err];

        tot_iter_number = tot_iter_number + 1;

        % Convergence Check
        if err < 1e-10        
         break;
       elseif isnan(err) || err > 100
          diverge = 1;
          break;  
        end
           errors(i) = err;

        %% Contraction rate calculation
           if i > 1
               contraction_rates(i - 1) = errors(i) / errors(i - 1);
           end
    end
    % Calculate contraction rates based on available data
    if length(contraction_rates) < 3
        arithmetic_mean = mean(contraction_rates);  % Use mean of all if less than 3 rates
    else
        arithmetic_mean = mean(contraction_rates(end - 2:end));  % Use mean of last two rates
    end
    %% Update current state
     S_current = S_iter;
     b_n = bs(S_current);
     iteration_no;
     EA;
end

iter_per_timestep = tot_iter_number / t_end;
