%% DD Toy-model - Error Vs Number of Iterations
% This script solves the DD Toy model using iterative 
% scheme (M scheme).

function [M_iter_per_step_values, h_M_values, iteration_M_iter, EA_M] = M_Error_MethodSolver(L, R, n_values, Mb, MB, T, dt, ee, m)

    %% Initialize arrays to store values
    M_iter_per_step_values = zeros(size(n_values));
    h_M_values = zeros(size(n_values));
    iteration_M_iter = zeros(size(n_values));

    for n_index = 1:length(n_values)

        n = n_values(n_index);
        dx = (R - L) / n;
        h = 1 / dx;
        h_M_values(n_index) = h;

        x = zeros(1, n + 1);
        for i = 1:n + 1
            x(i) = L + dx * (i - 1);
        end

        %% Number of Time Steps
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
        EA_M = [];
        diverge = 0; 
        iteration_no = 50;
        tot_iter_number = 0;

        %% Start of time loop
        for j = t_start:t_end 

            t(j) = (j - 1) * dt;

            if rem(j, 1) == 0
                EA_M = [];
            end

            S_iter = S_current;

            %% Precomputing arrays
            tm = (dt / dx^2) .* ones(1, n - 1);
            tp = (dt / dx^2) .* ones(1, n - 1);
            Main_D = [1, (tm + tp), 1];
            Lower_D = [0, -tm, 0];
            Upper_D = [0, -tp, 0];
           

            %% Iteration loop
            for i = 1:iteration_no

                % Compute b, B, db/ds, and dB/ds values
                b_iter = bs(S_iter);
                B_iter = BBs(S_iter);
                dBds_iter = dBBds(S_iter);
                dbds_iter = dbds(S_iter);

                % Define Lb and LB
                Lb = min(max(dbds_iter + Mb * dt, 2 * Mb * dt), 1);
                LB = min(max(dBds_iter + MB * dt, 2 * MB * dt), 1);

                % Matrix computation
                AA = Main_D + (Lb ./ LB) .* (1 - ee * dt);
                AA(1) = 1; 
                AA(end) = 1;

                f = b_n(2:n) - (1 - ee * dt) .* b_iter(2:n) + ...
                    (1 - ee * dt) .* B_iter(2:n) .* (Lb(2:n) ./ LB(2:n));

                F = [Pl, f, Pr];

                if i == 1
                    W_intiter = B_iter;
                else
                    W_intiter = W_iter;
                end

                W_iter = tridiagQQ(AA, Lower_D, Upper_D, F)';

                S_intiter = S_iter;
                
                S_iter = (1 ./ LB) .* (W_iter - B_iter) + S_intiter;

                dW_dx = gradient(W_iter - W_intiter, dx);

                %% Error calculation
                err = sqrt(trapz(x, (S_iter - S_intiter).^2 + dt * (dW_dx.^2)));

                EA_M = [EA_M, err];
               
                tot_iter_number = tot_iter_number + 1;

                %% Convergence check
                if err < 1e-10
                    break;
                elseif isnan(err) || err > 100
                    diverge=1;
                    fprintf('Divergence detected at MScheme %d\n', i);
                    break;
                end

            end 
           if diverge ==1
              iter_per_timestep = inf;
           end
            %% Variable Update
            S_current = S_iter;
            b_n = bs(S_current);
        end

        %% Compute iteration per timestep
        if diverge~= 1
        
        iter_per_timestep = tot_iter_number / t_end;
        end
        M_iter_per_step_values(n_index) = iter_per_timestep;
        iteration_M_iter(n_index) = i;
        h_M_values(n_index) = h;

    end
end
