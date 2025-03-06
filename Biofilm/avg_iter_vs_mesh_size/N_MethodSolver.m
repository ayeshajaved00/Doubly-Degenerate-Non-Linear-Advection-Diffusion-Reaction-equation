%% Biofilm Model - Average iteration Vs Error
% This script solves the Biofilm model using iterative 
% scheme (Newton scheme).

function [N_iter_per_step_values, h_N_values] = N_MethodSolver(L, R, n_values, T, dt, ee, m)

    % Initialize an array to store iter_per_step for each n
    N_iter_per_step_values = zeros(size(n_values));
    h_N_values = zeros(size(n_values));
    iteration_N_iter = zeros(size(n_values));
    M = 0;

    for n_index = 1:length(n_values)
        n = n_values(n_index);
        dx = (R - L) / n;
        h = 1 / dx;
        x = linspace(L, R, n + 1);

        %% The Number of Time Steps
        N = floor(T / dt) + 1;
        t = zeros(1, N);

        %% Definition of Initial and Boundary Condition
        Pl = 0;
        Pr = 0;
        S0 = InitialWW(x, m, 0, 1);
        b_n = S0;
        S_current = S0;

        %% Basic Iteration Information
        t_start = 1;
        t_end = N;
        diverge = 0;
        iteration_no = 50;
        tot_iter_number = 0;

        %% Start of Time Loop  
        for j = t_start:t_end
            t(j) = (j - 1) * dt;
            EA_N = [];
            S_iter = S_current;

            %% Precomputing the Arrays
            tm = (dt / dx^2) * ones(1, n - 1);
            tp = (dt / dx^2) * ones(1, n - 1);

            for i = 1:iteration_no
                b_iter = bs(S_iter);
                B_iter = BBs(S_iter);
                dbds_iter = dbds(S_iter);
                dBds_iter = dBBds(S_iter);

                %% Defining LB and Lb
                Lb = min(max(dbds_iter + M * dt, 2 * M * dt), 1);
                LB = min(max(dBds_iter + M * dt, 2 * M * dt), 1);

                Main_D = [1 (1 - ee .* (1 - b_n(2:n)) .* dt) .* (Lb(2:n)) + ((tm + tp) .* LB(2:n)) 1];
                Lower_D = [0 -tm .* LB(1:n-1) 0];
                Upper_D = [0 -tp .* LB(3:n+1) 0];

                f = b_n(2:n) - (1 - ee .* (1 - b_n(2:n)) .* dt) .* b_iter(2:n) + ...
                    (1 - ee .* (1 - b_n(2:n)) .* dt) .* Lb(2:n) .* S_iter(2:n) + ...
                    ((dt / dx) .* (S_iter(2:n)) .* LB(2:n) .* (2 / dx) - ...
                    (dt / dx) .* (2 / dx) .* B_iter(2:n) - (dt / dx^2) .* LB(1:n-1) .* S_iter(1:n-1) + ...
                    B_iter(1:n-1) .* (dt / dx^2) - (dt / dx^2) .* LB(3:n+1) .* (S_iter(3:n+1)) + ...
                    (dt / dx^2) .* B_iter(3:n+1));

                F = [Pl f Pr];
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
                                
                err = sqrt(trapz(x,  Lb .* LB .*(S_iter - S_intiter).^2 + dt * (dW_dx.^2)));

                EA_N = [EA_N, err];
                
                tot_iter_number = tot_iter_number + 1;

                %% Convergence Check
                if err < 1e-6
                    break;
                elseif isnan(err) || err > 100
                    diverge = 1;
                    break;
                end
            end

            %% Variable Update
            S_current = S_iter;
            b_n = bs(S_current);

            if diverge == 1
                iter_per_timestep = inf;
                break;
            end
        end

       if diverge ~= 1
            iter_per_timestep = tot_iter_number / t_end;
       end
        N_iter_per_step_values(n_index) = iter_per_timestep;
        h_N_values(n_index) = h;
        iteration_N_iter(n_index) = i;
    end
end
