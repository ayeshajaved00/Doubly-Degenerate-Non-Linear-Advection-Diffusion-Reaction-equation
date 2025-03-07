%% Biofilm Model - Convergence Order
% This script solves the Biofilm model using iterative 
% scheme (N scheme) and calculates the convergence order of the scheme.

function [h_N_values, average_orders_N, orders, EA_N, errors] = N_Order_MethodSolver(L, R, n_values, T, dt, ee, m)

    % Initialize arrays to store results
    h_N_values = zeros(size(n_values));
    average_orders_N = zeros(size(n_values));
    M = 0;

    for n_index = 1:length(n_values)
        n = n_values(n_index);
        dx = (R - L) / n;
        h = 1 / dx;

        x = zeros(1, n + 1);
        for i = 1:n + 1
            x(i) = L + dx * (i - 1);
        end

        %% The Number of Time Steps
        N = floor(T / dt) + 1;
        t = zeros(1, N);

        %% Definition of initial and boundary conditions
        Pl = 0;
        Pr = 0;
        S0 = InitialWW(x, m, 0, 1);
        b_n = S0;
        S_current = S0;

        %% Basic iteration information
        t_start = 1;
        t_end = N;
        orders = [];
        errors = [];
        EA_N = [];
        diverge = 0;
        iteration_no = 50;
        tot_iter_number = 0;

        for j = t_start:t_end
            t(j) = (j - 1) * dt;

            if rem(j, 10) == 0
                EA_N = [];
                orders = [];
                errors = [];
            end

            S_iter = S_current;
            tm = (dt / dx^2) .* ones(1, n - 1);
            tp = (dt / dx^2) .* ones(1, n - 1);

            for i = 1:iteration_no
                b_iter = bs(S_iter);
                B_iter = BBs(S_iter);
                dbds_iter = dbds(S_iter);
                dBds_iter = dBBds(S_iter);

                Lb = min(max(dbds_iter + M * dt, 2 * M * dt), 1);
                LB = min(max(dBds_iter + M * dt, 2 * M * dt), 1);

                Main_D = [1, (1 - ee * (1 - b_n(2:n)) * dt) .* Lb(2:n) + ((tm + tp) .* LB(2:n)), 1];
                Lower_D = [0, -tm .* LB(1:n-1), 0];
                Upper_D = [0, -tp .* LB(3:n+1), 0];

                A(1) = 1;
                A(end) = 1;

                f = b_n(2:n) - (1 - ee .* (1 - b_n(2:n)) * dt) .* b_iter(2:n) + ...
                    (1 - ee .* (1 - b_n(2:n)) * dt) .* Lb(2:n) .* S_iter(2:n) + ...
                    ((dt / dx) .* (S_iter(2:n)) .* LB(2:n) .* (2 / dx) - ...
                    (dt / dx) .* (2 / dx) .* B_iter(2:n) - (dt / dx^2) .* LB(1:n-1) .* S_iter(1:n-1) + ...
                    B_iter(1:n-1) .* (dt / dx^2) - (dt / dx^2) .* LB(3:n+1) .* (S_iter(3:n+1)) + ...
                    (dt / dx^2) .* B_iter(3:n+1));

                F = [Pl, f, Pr];
                S_intiter = S_iter;
                S_iter = tridiagQQ(Main_D, Lower_D, Upper_D, F)';

                if i == 1
                    W_intiter = B_iter;
                else
                    W_intiter = W_iter;
                end
                W_iter = LB .* (S_iter - S_intiter) + B_iter;

                %% Error calculation
                dW_dx = gradient(W_iter - W_intiter, dx);
                
                err = sqrt(trapz(x,  (S_iter - S_intiter).^2  + dt * (dW_dx.^2)));

                EA_N = [EA_N, err];
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

                %% Convergence order
                for k = 3:length(errors)
                    orders(k - 2) = abs(log(errors(k) / errors(k - 1)) / log(errors(k - 1) / errors(k - 2)));
                end
            end

             if diverge == 1
                    Convergence_order = inf; 
             end
             if diverge ~= 1
                 if length(orders) < 3
                    Convergence_order = mean(orders);  % Use mean of all if less than 3 rates
                 else
                    Convergence_order = mean(orders(end - 2:end));  % Use mean of last two rates
                 end
             end
            %% Variable Update
            S_current = S_iter;
            b_n = bs(S_current);
        end

         h_N_values(n_index) = h;
        average_orders_N(n_index) = Convergence_order;
    end
end
