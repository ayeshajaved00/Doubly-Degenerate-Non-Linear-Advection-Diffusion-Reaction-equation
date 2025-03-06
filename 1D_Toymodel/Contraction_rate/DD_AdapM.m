%% DD Toymodel - Contraction rate
% This script solves the Richards equation 
% using an iterative scheme (M-Adaptive scheme) 
% and calculates the contraction rate.

clc;
clear all;
close all;

%% Input Values
L = input('Left starting point of the domain: ');
R = input('Right starting point of the domain: ');
n = input('Number of cells: ');
T = input('Total runtime: ');
dt = input('Time step: ');

%% Adaptivity Parameters
adap_num_values = 25;  % Number of values in the range
adap_start_value = -10; % Starting value
adap_end_value = -2;    % Ending value
adap_step_size = (adap_end_value - adap_start_value) / (adap_num_values - 1);

ee = 0.5; % Source term
m = 6;
MB = 1;
Mb = 1;

dx = (R - L) / n;

x = zeros(1,n+1);

for i = 1:n+1
   x(i) = L + dx * (i-1); 
end

%% The Number of Time Steps
N = floor(T/dt) + 1;  
t = zeros(1, N); 

%% Definition of initial and boundary conditions
Pl = 0; 
Pr = 0;

% Initial condition
S0 = InitialWW(x, m, 0, 1);
b_n = S0;
S_current = S0;

%% Basic iteration information
t_start = 1; 
t_end = N;
diverge = 0;
iteration_no = 50;    
EA = [];
contraction_rates = [];
errors = [];

%% Start of time loop
tot_iter_number = 0;

for j = t_start:t_end
    t(j) = (j-1) * dt;

    if rem(j, 10) == 0
        EA = [];
        contraction_rates = [];
        errors = [];
    end

    S_iter = S_current;

    %% Precomputing the arrays
    tm = ((dt/dx^2)) * ones(1, n-1);
    tp = ((dt/dx^2)) * ones(1, n-1);

    Main_D = [1, (tm + tp), 1];  
    Lower_D = [0, -tm, 0];
    Upper_D = [0, -tp, 0];

    Mb(1) = Mb(end); 
    Mb = [Mb(1)];
    MB(1) = MB(end); 
    MB = [MB(1)];

    % Initial variable updates
    b_iter = bs(S_iter); % This is the value of b at S^{i-1}_n
    B_iter = BBs(S_iter); % This is the value of B at S^{i-1}_n
    dBds_iter = dBBds(S_iter); % This is the value of B' at S^{i-1}_n
    dbds_iter = dbds(S_iter); % This is the value of b' at S^{i-1}_n

    for i = 1:iteration_no      
        %% Definition of L_b and L_B
        Lb = min(max(dbds_iter + Mb(i) * dt, 2 * Mb(i) * dt), 1);
        LB = min(max(dBds_iter + MB(i) * dt, 2 * MB(i) * dt), 1);

        %% Matrix computation
        AA = Main_D + (Lb ./ LB) .* (1 - ee * dt);
        AA(1) = 1; 
        AA(end) = 1;

        f = b_n(2:n) - (1 - ee * dt) .* b_iter(2:n) + (1 - ee * dt) .* B_iter(2:n) .* ((Lb(2:n)) ./ (LB(2:n)));
        F = [Pl, f, Pr];

        if i == 1
            W_intiter = B_iter;
        else
            W_intiter = W_iter;
        end

        W_iter = tridiagQQ(AA, Lower_D, Upper_D, F)';
        S_intiter = S_iter;
        S_iter = (1 ./ LB) .* (W_iter - B_iter) + S_intiter;

        U_iter = b_iter + Lb .* (S_iter - S_intiter);
        b_iter = bs(S_iter); % This is the value of b at S^{i-1}_n
        B_iter = BBs(S_iter); % This is the value of B at S^{i-1}_n
        dBds_iter = dBBds(S_iter); % This is the value of B' at S^{i-1}_n
        dbds_iter = dbds(S_iter); % This is the value of b' at S^{i-1}_n

        %% Error Calculation
        dW_dx = gradient(W_iter - W_intiter, dx);
       
        err = sqrt(trapz(x, (S_iter - S_intiter).^2 + dt * (dW_dx.^2)));

        EA = [EA, err];
        tot_iter_number = tot_iter_number + 1;

        %% Convergence Check
        if err < 1e-10
            break;
        elseif isnan(err) || err > 100
            diverge = 1;
            break;
        end
        errors(i) = err;

        %% Adaptivity
        if rem(i, 3) == 2
            U = U_iter - b_iter;
            O = W_iter - B_iter;

            for j1 = 1:adap_num_values
                value = adap_start_value + (j1 - 1) * adap_step_size;
                value = 10^value;
                Lb = min(max(dbds_iter + value * dt, 2 * value * dt), 1);
                LB = min(max(dBds_iter + value * dt, 2 * value * dt), 1);
                G = eta(LB ./ Lb, U, O, dx);
                if G < err
                    break;
                end
            end

            Mb(i + 1) = value;
            MB(i + 1) = value;
        else
            Mb(i + 1) = Mb(i);
            MB(i + 1) = MB(i);
        end

        %% Contraction Rate Calculation
        if i > 1
            contraction_rates(i - 1) = errors(i) / errors(i - 1);
        end
    end

    %% Variable Update
    S_current = S_iter;
    b_n = bs(S_current);
end

%% Post-Processing
if length(contraction_rates) < 3
    arithmetic = mean(contraction_rates);
else
    arithmetic = mean(contraction_rates(end - 2:end));
end
iter_per_timestep = tot_iter_number / t_end;
