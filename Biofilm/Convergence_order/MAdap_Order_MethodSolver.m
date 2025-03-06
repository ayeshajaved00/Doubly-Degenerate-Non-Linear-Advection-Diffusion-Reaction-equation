%% Biofilm Model - Convergence Order
% This script solves the Biofilm model using iterative 
% scheme (M-Adaptive scheme) and calculate the Convergence order.

function [h_MAdap_values, average_orders_AdapM, EA_MAdap, orders] = MAdap_Order_MethodSolver(L, R, n_values, T, dt, ee, m)

% Initialize an array to store iter_per_step for each n
h_MAdap_values = zeros(size(n_values));
average_orders_AdapM = zeros(size(n_values));
Mb = 1; MB = 1;

% Adaptivity parameters
adap_num_values = 25; % Number of values in the range
adap_start_value = -10; % Starting value
adap_end_value = -2; % Ending value
adap_step_size = (adap_end_value - adap_start_value) / (adap_num_values - 1);

for n_index = 1:length(n_values)

    n = n_values(n_index);
    dx = (R - L) / n;
    h = 1 / dx;
    x = zeros(1, n + 1);

    for i = 1:n + 1
        x(i) = L + dx * (i - 1);
    end

    % The Number of Time Steps
    N = floor(T / dt) + 1;
    t = zeros(1, N);

    % Boundary Condition
    Pl = 0;
    Pr = 0;

    % Initial condition
    S0 = InitialWW(x, m, 0, 1);
    b_n = S0;
    S_current = S0;

    % Basic iteration information
    t_start = 1; t_end = N;
    EA_MAdap = []; errors = []; orders = []; Convergence_order=[];
    diverge = 0;
    iteration_no = 50;

    % Start of time loop
    tot_iter_number = 0;

      for j = t_start:t_end

        t(j) = (j - 1) * dt;

        if rem(j, 10) == 0
            EA_MAdap = []; errors = []; orders = [];
        end

        S_iter = S_current;

        % Precomputing the arrays
        tm = ((dt / dx^2)) .* ones(1, n - 1);
        tp = ((dt / dx^2)) .* ones(1, n - 1);

        Main_D = [1 (tm + tp) 1];  
        Lower_D = [0 -tm 0];
        Upper_D = [0 -tp 0];

        Mb(1) = Mb(end); Mb = [Mb(1)];
        MB(1) = MB(end); MB = [MB(1)];

        % Basic variable update
        b_iter = bs(S_iter); % This is the value of b at S^{i-1}_n
        B_iter = BBs(S_iter); % This is the value of B at S^{i-1}_n
        dbds_iter = dbds(S_iter); % This is the value of b' at S^{i-1}_n
        dBds_iter = dBBds(S_iter); % This is the value of B' at S^{i-1}_n

           for i = 1:iteration_no

            % Defining LB and Lb

            % M-scheme
            Lb = min(max(dbds_iter + Mb(i) * dt, 2 * Mb(i) * dt), 1);
            LB = min(max(dBds_iter + MB(i) * dt, 2 * MB(i) * dt), 1);

            % AA Matrix and forcing term
            AA = Main_D + (Lb ./ LB) .* (1 - ee * (1 - b_n) * dt);

            f = b_n(2:n) - (1 - ee .* (1 - b_n(2:n)) * dt) .* b_iter(2:n) +...
                (1 - ee .* (1 - b_n(2:n)) * dt) .* B_iter(2:n) .* ((Lb(2:n)) ./ (LB(2:n)));

            F = [Pl f Pr];
            if i == 1
                W_intiter = B_iter;
            else
                 W_intiter = W_iter;
            end
            %% Solve the linear system
            W_iter = tridiagQQ(AA, Lower_D, Upper_D, F)';

            S_intiter = S_iter;
            S_iter = (1 ./ LB) .* (W_iter - B_iter) + S_intiter;
            U_iter = b_iter + Lb .* (S_iter - S_intiter);

            %% Update function
            b_iter = bs(S_iter);
            B_iter = BBs(S_iter);
            dbds_iter = dbds(S_iter);
            dBds_iter = dBBds(S_iter);

            %% Error calculation
            dW_dx = gradient(W_iter - W_intiter, dx);

            err = sqrt(trapz(x,  (S_iter - S_intiter).^2  + dt *(dW_dx.^2)));

            EA_MAdap = [EA_MAdap, err];    

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

               
            %% Adaptivity
            if rem(i, 2) == 1
                U = U_iter - b_iter;
                O = W_iter - B_iter;

                for j1 = 0:adap_num_values
                    value = adap_start_value + (j1) * adap_step_size;
                    value = 10.^value;
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
         EA_MAdap;
  
      end
    h_MAdap_values(n_index) = h;
    average_orders_AdapM(n_index) = Convergence_order;
end
end
