%% DD Toy-model - Convergence Order
% This script solves the DD Toy model using iterative 
% schemes (M scheme) and calculates the convergence order of the scheme.

function [h_M_values, average_orders_M] = M_Order_MethodSolver(L, R, n_values, Mb, MB, T, dt, ee, m)
    %% Initialize arrays to store h values and average orders for each n
    h_M_values = zeros(size(n_values));
    average_orders_M = zeros(size(n_values));

    for n_index = 1:length(n_values)
        n = n_values(n_index);
        dx = (R - L) / n;
        h = 1 / dx;

        %% Initialize spatial grid
        x = zeros(1, n + 1);
        for i = 1:n + 1
            x(i) = L + dx * (i - 1);
        end

        %% Number of time steps
        N = floor(T / dt) + 1;
        t = zeros(1, N);

        %% Define initial and boundary conditions
        Pl = 0; Pr = 0;
        S0 = InitialWW(x, m, 0, 1); % Initial condition
        b_n = S0;
        S_current = S0;

        %% Basic iteration information
        iteration_no = 50;
        EA = [];
        orders = [];
        errors = [];
        diverge = 0;

        tot_iter_number = 0;

        %% Time loop
        for j = 1:N
            t(j) = (j - 1) * dt;

            if rem(j, 1) == 0
                EA = [];
                orders = [];
                errors = [];
            end

            S_iter = S_current;

            %% Precomputing arrays
            tm = (dt / dx^2) * ones(1, n - 1);
            tp = (dt / dx^2) * ones(1, n - 1);

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
                AA = Main_D + (Lb ./ LB).* (1 - ee * dt);
                AA(1) = 1; AA(end) = 1;

                f = b_n(2:n) - (1 - ee * dt).* b_iter(2:n) + ...
                    (1 - ee * dt).* B_iter(2:n).* (Lb(2:n)./ LB(2:n));
                F = [Pl, f, Pr];
                if i == 1
                   W_intiter = B_iter;
               else
                   W_intiter = W_iter;
               end
                W_iter = tridiagQQ(AA, Lower_D, Upper_D, F)';

                S_intiter = S_iter;
                S_iter = (1 ./ LB).* (W_iter - B_iter) + S_intiter;

                %% Error calculation
                dW_dx = gradient(W_iter - W_intiter, dx);
                err = sqrt(trapz(x,(S_iter - S_intiter).^2 + dt * (dW_dx.^2)));      
                EA = [EA, err];

                tot_iter_number = tot_iter_number + 1;
                errors(i) = err;

                %% Convergence check
                if err < 1e-10
                    break;
                elseif isnan(err) || err > 100
                    diverge = 1;
                    break;
                end

                %% Compute convergence orders
                for k = 3:length(errors)
                    orders(k - 2) = abs(log(errors(k) / errors(k - 1)) / log(errors(k - 1) / errors(k - 2)));
                end
            end
            if diverge == 1
                Convergence_order =inf;
            end
            %% Update variables
            S_current = S_iter;
            b_n = bs(S_current);
        end

        % Compute average convergence order
        if diverge ~=1
        if length(orders) < 3
            Convergence_order = mean(orders);
        else
            Convergence_order = mean(orders(end - 2:end));
        end
        end

        h_M_values(n_index) = h;
        average_orders_M(n_index) = Convergence_order;
    end
end
