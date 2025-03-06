function [h_M_values, average_orders_M] = M_Order_MethodSolver(L, R, n_values, b_values, B_values, Mbb, MBB, T, dt, ee, m)

    global Ks
    global lambda

    %% Initialize arrays to store results
    h_M_values = zeros(size(n_values));
    average_orders_M = zeros(size(n_values));

    for n_index = 1:length(n_values)
        n = n_values(n_index);
        dx = (R - L) / n;
        h = 1 / dx;
        h_M_values(n_index) = h;

        %% Initialize spatial grid
            x = zeros(1,n+1);

          for i = 1:n+1
             x(i) = L+dx*(i-1); 
          end
        %% Number of time steps
        N = floor(T / dt) + 1;
        t = zeros(1, N);

        %% Boundary conditions
        Pl = 0; Pr = 0;

        %% Initial condition
        S0 = InitialWW(x, m, 0, 1);
        b_n = S0;
        S_current = S0;

        %% Basic iteration information
        t_start = 1; t_end = N;
        EA_M = []; orders = []; errors = [];
        diverge = 0; iteration_no = 50;

        %% Start of time loop
        tot_iter_number = 0;

        for j = t_start:t_end
            t(j) = (j - 1) * dt;

            if rem(j, 10) == 0
                EA_M = []; orders = []; errors = [];
            end

            S_iter = S_current;

            %% Precompute the arrays
            tm = (dt / dx^2) * ones(1, n - 1);
            tp = (dt / dx^2) * ones(1, n - 1);
            Main_D = [1, (tm + tp), 1];
            Lower_D = [0, -tm, 0];
            Upper_D = [0, -tp, 0];

            for i = 1:iteration_no
                S_iter = max(S_iter, 0);
                index = idx(S_iter, b_values);
                b_iter = bs(S_iter, index, b_values);
                B_iter = BB(S_iter, index, B_values);
                dbds_iter = dbds(S_iter, index, b_values);
                dBds_iter = dBBds(S_iter, index, B_values);

                %% Advection term
                F_iter = Ks(b_iter, lambda);

                %% Define Lb and LB (M-scheme)
                Lb = min(max(dbds_iter + Mbb * dt, 2 * Mbb * dt), 1);
                LB = min(max(dBds_iter + MBB * dt, 2 * MBB * dt), 1);

                %% Assemble system matrix and RHS
                AA = Main_D + (Lb ./ LB) .* (1 - ee * dt);
                AA(1) = 1; AA(end) = 1;

                f = b_n(2:n) - (1 - ee * dt) .* b_iter(2:n) + ...
                    (1 - ee * dt) .* B_iter(2:n) .* (Lb(2:n) ./ LB(2:n)) + ...
                    (dt / dx) * (F_iter(1:n - 1) - F_iter(2:n));
                F = [Pl, f, Pr];

                %% Solve the tridiagonal system
                if i == 1
                   W_intiter = B_iter;
                else
                   W_intiter = W_iter;
                end
                W_iter = tridiagQQ(AA, Lower_D, Upper_D, F)';

                S_intiter = S_iter;
                S_iter = (1 ./ LB) .* (W_iter - B_iter) + S_intiter;

                %% Error calculation
                dW_dx = gradient(W_iter - W_intiter, dx);
                
                err = sqrt(trapz(x,  (S_iter - S_intiter).^2  + dt * (dW_dx.^2)));

                EA_M = [EA_M, err];
                
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

                %% Order of convergence
                for k = 3:length(errors)
                    orders(k - 2) = abs(log(errors(k) / errors(k - 1)) / log(errors(k - 1) / errors(k - 2)));
                end
            end

            %% Update variables
            S_current = S_iter;
            b_n = bs(S_current, index, b_values);
            EA_M;
            if diverge == 1
                Convergence_order = inf;
                break;
            end
        end

        %% Final convergence order calculation
        if diverge ~= 1
            Convergence_order = mean(orders(end - 2:end));
        end

        h_M_values(n_index) = h;
        average_orders_M(n_index) = Convergence_order;
    end
end
