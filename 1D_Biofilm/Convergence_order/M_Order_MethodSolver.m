%% Biofilm Model - Convergence Order
% This script solves the Biofilm model using iterative 
% scheme (M scheme) and calculate the Convergence order.

function [h_M_values,average_orders_M,orders] = M_Order_MethodSolver(L, R, n_values, Mbb, MBB, T, dt, ee, m)

    %% Initialize arrays to store results for each n
    h_M_values = zeros(size(n_values));
    average_orders_M = zeros(size(n_values));

    for n_index = 1:length(n_values)
        n = n_values(n_index);
        
        % Calculate step size
        dx = (R - L) / n;
        
        % Calculate h
        h = 1 / dx;

        % Initialize x vector
        x = zeros(1, n + 1);
        
        % Generate x values
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

        EA_M = [];
        errors = [];
        orders = [];
        iteration_no = 50;
        diverge = 0;

        %% Start of time loop
        tot_iter_number = 0;

        for j = t_start:t_end
            t(j) = (j - 1) * dt;

            if rem(j, 10) == 0
                EA_M = [];   
                errors = []; 
                orders = []; 
            end

            S_iter = S_current;

            %% Precomputing the arrays
            tm = (dt / dx^2) .* ones(1, n - 1);
            tp = (dt / dx^2) .* ones(1, n - 1);
            Main_D = [1, (tm + tp), 1];  
            Lower_D = [0, -tm, 0];
            Upper_D = [0, -tp, 0];

            for i = 1:iteration_no
                b_iter = bs(S_iter); % This is the value of b at S^{i-1}_n
                B_iter = BBs(S_iter); % This is the value of B at S^{i-1}_n
                dbds_iter = dbds(S_iter); % This is the value of b' at S^{i-1}_n
                dBds_iter = dBBds(S_iter); % This is the value of B' at S^{i-1}_n

                %% Defining LB and Lb
                % M-scheme
                Lb = min(max(dbds_iter + Mbb * dt, 2 * Mbb * dt), 1);
                LB = min(max(dBds_iter + MBB * dt, 2 * MBB * dt), 1);

                AA = Main_D + (Lb ./ LB) .* (1 - ee * (1 - b_n) * dt);

                f = b_n(2:n) - (1 - ee .* (1 - b_n(2:n)) * dt) .* b_iter(2:n) +...
                    (1 - ee .* (1 - b_n(2:n)) * dt) .* B_iter(2:n) .* ((Lb(2:n)) ./ (LB(2:n)));

                F = [Pl, f, Pr];
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
               
                err = sqrt(trapz(x,  (S_iter - S_intiter).^2  + dt *(dW_dx.^2)));

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
                    Convergence_order = mean(orders(end - 2:end)); 
                 end
             end
            %% Variable Update
            S_current = S_iter;
            b_n = bs(S_current);
            EA_M;

            if diverge == 1
                Convergence_order = inf;
                break;
            end
             if diverge ~= 1
                Convergence_order = mean(orders(end - 2:end));
             end
        end

        h_M_values(n_index) = h;
        average_orders_M(n_index) = Convergence_order;
    end

end
