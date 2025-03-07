clc;
clear all;
close all;

%% Input Values
L = input('Left starting point of the domain: ');
R = input('Right starting point of the domain: ');
n = input('Number of cells: ');
T = input('Total runtime: ');
dt = input('Time step: ');

m = 6; ee = 0;
M = 0;
dx = (R - L) / n;

x = zeros(1, n + 1);
for k = 1:n + 1
    x(k) = L + dx * (k - 1);
end

%% Number of Time Steps
N = floor(T / dt) + 1;
t = zeros(1, N);

%% Initial and Boundary Conditions
Pl = 0; Pr = 0; % Boundary conditions
S0 = InitialWW(x, m, 0, 1); % Initial condition

b_n = S0;
S_current = S0;

%% Basic Iteration Information
t_start = 1; t_end = N;
EA = []; contraction_rates = []; errors = []; orders = [];
iteration_no = 50; diverge = 0;

%% Time Loop
tot_iter_number = 0;

for j = t_start:t_end
    no_it_pertimestep = 0;
    t(j) = (j - 1) * dt;

    if rem(j, 10) == 0
        EA = []; contraction_rates = []; errors = []; orders = [];
    end

    S_iter = S_current;

    %% Precomputing Arrays
    tm = ((dt / dx^2)) .* ones(1, n - 1);
    tp = ((dt / dx^2)) .* ones(1, n - 1);

    for i = 1:iteration_no
        b_iter = bs(S_iter, m); % b(S^{i-1}_n)
        B_iter = BBs(S_iter, m); % B(S^{i-1}_n)
        dBds_iter = dBBds(S_iter, m); % B'(S^{i-1}_n)
        dbds_iter = dbds(S_iter, m); % b'(S^{i-1}_n)

        %% Definition of Lb and LB
        Lb = min(max(dbds_iter + M * dt, 2 * M * dt), 1);
        LB = min(max(dBds_iter + M * dt, 2 * M * dt), 1);

        Main_D = [1, (1 - ee * dt) .* Lb(2:n) + ((tm + tp) .* LB(2:n)), 1];
        Lower_D = [0, -tm .* LB(1:n - 1), 0];
        Upper_D = [0, -tp .* LB(3:n + 1), 0];

        Main_D(1) = 1; Main_D(end) = 1;

        f = b_n(2:n) - (1 - ee * dt) .* b_iter(2:n) ...
            + (1 - ee * dt) .* Lb(2:n) .* S_iter(2:n) ...
            + ((dt / dx) .* (S_iter(2:n)) .* LB(2:n) .* (2 / dx)) ...
            - ((dt / dx) * (2 / dx) .* B_iter(2:n)) ...
            - ((dt / dx^2) .* LB(1:n - 1) .* S_iter(1:n - 1)) ...
            + B_iter(1:n - 1) .* (dt / dx^2) ...
            - ((dt / dx^2) .* LB(3:n + 1) .* S_iter(3:n + 1)) ...
            + ((dt / dx^2) .* B_iter(3:n + 1));

        F = [(1 + (LB(1) ./ Lb(1))) * Pl, f, (1 + (LB(end) ./ Lb(end))) * Pr];

        S_intiter = S_iter;
        S_iter = tridiagQQ(Main_D, Lower_D, Upper_D, F)';

        if i == 1
            W_intiter = B_iter;
        else
            W_intiter = W_iter;
        end
        W_iter = LB .* (S_iter - S_intiter) + B_iter;

        %% Error Calculation
         dW_dx = gradient(W_iter - W_intiter, dx);
        
         err = sqrt(trapz(x,(S_iter - S_intiter).^2  + dt *(dW_dx.^2)));
        
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

        %% Contraction Rate
        if i > 1
            contraction_rates(i - 1) = errors(i) / errors(i - 1);
        end
    end
  %% Stop Simulation if Diverged
    if diverge == 1
        fprintf('Simulation stopped due to divergence at time step %d.\n', j);
        break; % Exit the outer loop
    end
    %% Update Variables
    S_current = S_iter;
    b_n = bs(S_current, m);


end
  %% Calculate Arithmetic Mean Contraction Rate
    if length(contraction_rates) < 3
        arithmetic = mean(contraction_rates);
    else
        arithmetic = mean(contraction_rates(end - 2:end));
    end
iter_per_timestep = tot_iter_number / t_end;
