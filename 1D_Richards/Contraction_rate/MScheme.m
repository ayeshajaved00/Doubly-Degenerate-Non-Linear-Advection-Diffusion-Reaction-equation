%% Richards Equation - Contraction rate
% This script solves the Richards equation using iterative 
% scheme (M-scheme) and visualizes the contraction rate

clc;
clear all;
close all;

Richards;

%% Input values
L = input('Left starting point of the domain: ');
R = input('Right starting point of the domain: ');
n = input('Number of cells: ');
T = input('Total runtime: ');
dt = input('Time step: ');

dx = (R - L) / n;

Mb = 0.01;
MB = 0.01;

m = 6; ee = 0.1; 

x = zeros(1, n + 1);
for l = 1:n + 1
    x(l) = L + dx * (l - 1);
end

%% The Number of Time Steps
N = floor(T / dt) + 1;
t = zeros(1, N);

%% Boundary Condition
Pl = 0;
Pr = 0;

%% Initial condition
S0 = InitialWW(x, m, 0, 1);
b_n = S0;
S_current = S0;

%% Basic iteration information
t_start = 1; 
t_end = N;

EA = [];  
contraction_rates = []; 
errors = []; 
diverge = 0;

iteration_no = 50; 
Array_iterations = [];  

%% Start of time loop
tot_iter_number = 0;
tic
for j = t_start:t_end 
    no_it_pertimestep = 0;
    t(j) = (j - 1) * dt;
    
    if rem(j, 10) == 0 
        EA = []; 
         errors = []; 
         contraction_rates = []; 
    end

    S_iter = S_current;

    %% Precomputing the arrays
    tm = (dt / dx^2) * ones(1, n - 1);
    tp = (dt / dx^2) * ones(1, n - 1);

    Main_D = [1, (tm + tp), 1];
    Lower_D = [0, -tm, 0];
    Upper_D = [0, -tp, 0];

    for i = 1:iteration_no 
        S_iter = max(S_iter, 0);

        index = idx(S_iter, b_values);
        b_iter = bs(S_iter, index, b_values); % Value of b at S^{i-1}_n
        B_iter = BB(S_iter, index, B_values); % Value of B at S^{i-1}_n
        dbds_iter = dbds(S_iter, index, b_values); % Value of b' at S^{i-1}_n
        dBds_iter = dBBds(S_iter, index, B_values); % Value of B' at S^{i-1}_n
        F_iter = Ks(b_iter, lambda); % Advection term

        %% Defining LB and Lb
        % M-scheme
        Lb = min(max(dbds_iter + Mb * dt, 2 * Mb * dt), 1);
        LB = min(max(dBds_iter + MB * dt, 2 * MB * dt), 1);

        AA = Main_D + (Lb ./ LB) .* (1 - ee * dt);
        AA(1) = 1; 
        AA(end) = 1;

        f = b_n(2:n) - (1 - ee * dt) .* b_iter(2:n) + ...
            (1 - ee * dt) .* B_iter(2:n) .* ((Lb(2:n)) ./ (LB(2:n))) + ...
            (dt / dx) * (F_iter(1:n-1) - F_iter(2:n));

        F = [Pl, f, Pr];
        if i == 1
           W_intiter = B_iter;
        else
           W_intiter = W_iter;
        end
        W_iter = tridiagQQ(AA, Lower_D, Upper_D, F)';
        S_intiter = S_iter;
        S_iter = (1 ./ LB) .* (W_iter - B_iter) + S_intiter;

        %% Error check
         dW_dx = gradient(W_iter - W_intiter, dx);
         
         err = sqrt(trapz(x,  (S_iter - S_intiter).^2 + dt *  (dW_dx.^2)));

         EA = [EA, err];   
        
         tot_iter_number = tot_iter_number + 1;

        %% Convergence Check
        if err < 1e-10
            break;
        elseif isnan(err) || err > 100
            diverge = 1;
            fprintf('Iteration %d: err = %e\n', i, err);
            fprintf('The solution diverged at time step %d, iteration %d.\n', j, i);
            break; 
        end

        errors(i) = err;

        if i > 1
            contraction_rates(i - 1) = (errors(i) / errors(i-1));
        end
    end
     if length(contraction_rates) < 3
        arithmetic = mean(contraction_rates);  % Use mean of all if less than 3 rates
    else
        arithmetic = mean(contraction_rates(end - 2:end)); 
     end

    %% Variable Update
    S_current = S_iter;
    b_n = bs(S_current, index, b_values);
    iteration_no;
    EA;
end
toc
iter_per_timestep = tot_iter_number / t_end;
