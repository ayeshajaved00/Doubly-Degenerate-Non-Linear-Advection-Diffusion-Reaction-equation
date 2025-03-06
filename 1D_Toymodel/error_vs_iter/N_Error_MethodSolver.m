%% DD Toy-model - Error Vs Number of Iterations
% This script solves the DD Toy model using iterative 
% scheme (Newton scheme).

function [N_iter_per_step_values, h_N_values, iteration_N_iter, EA_N] = N_Error_MethodSolver(L, R, n_values, T, dt, ee, m)
    M = 0;
    
    N_iter_per_step_values = zeros(size(n_values));
    h_N_values = zeros(size(n_values));
    iteration_N_iter = zeros(size(n_values));

    for n_index = 1:length(n_values)
        n = n_values(n_index);
        dx = (R - L) / n;
        h = 1 / dx;
        
        x = zeros(1, n + 1);
        for k = 1:n + 1
            x(k) = L + dx * (k - 1);
        end

        %% The Number of Time Steps 
        N = floor(T / dt) + 1;  
        t = zeros(1, N); 

        %% Definition of initial and boundary condition
        % Boundary Condition
        Pl = 0; 
        Pr = 0;

        % Initial condition
        S0 = InitialWW(x, m, 0, 1);
        S_current = S0;
        b_n = S0; 

        %% Basic iteration information
        t_start = 1; 
        t_end = N;
        diverge = 0;
        iteration_no = 50;     
        EA_N = [];

        %% Start of time loop  
        tot_iter_number = 0;

        for j = t_start:t_end 
            t(j) = (j - 1) * dt;
            
            if rem(j, 1) == 0 
                EA_N = [];
            end

            S_iter = S_current;

           %% Precomputing Arrays
          tm = (dt / dx^2) .* ones(1, n - 1);
          tp = (dt / dx^2) .* ones(1, n - 1);

         for i = 1:iteration_no 
            b_iter = bs(S_iter);        % b(S^{i-1}_n)
            B_iter = BBs(S_iter);       % B(S^{i-1}_n)
            dBds_iter = dBBds(S_iter); % B'(S^{i-1}_n)
            dbds_iter = dbds(S_iter);  % b'(S^{i-1}_n)

            %% Define Lb and LB
            Lb = min(max(dbds_iter + M * dt, 2 * M * dt), 1);
            LB = min(max(dBds_iter + M * dt, 2 * M * dt), 1);

            Main_D = [1, (1 - ee * dt) .* (Lb(2:n)) + ((tm + tp) .* LB(2:n)), 1];
            Lower_D = [0, -tm .* LB(1:n - 1), 0];
            Upper_D = [0, -tp .* LB(3:n + 1), 0];
             
            Main_D(1) = 1;  Main_D(end) = 1;

            f = b_n(2:n) - (1 - ee * dt).* b_iter(2:n) + ...
                (1 - ee * dt) .* Lb(2:n).* S_iter(2:n) + ...
                ((dt / dx) .* (S_iter(2:n)).* LB(2:n) .* (2 / dx) - ...
                (dt / dx) .* (2 ./ dx) .* B_iter(2:n) - ...
                (dt / dx^2) .* LB(1:n - 1) .* S_iter(1:n - 1) + ...
                B_iter(1:n - 1) .* (dt ./ dx^2) - ...
                (dt / dx^2) .* LB(3:n + 1) .* (S_iter(3:n + 1)) + ...
                (dt / dx^2) .* B_iter(3:n + 1));

            F = [(1 + (LB(1) / Lb(1))) * Pl, f, (1 + (LB(end) / Lb(end))) * Pr];

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
                
                err = sqrt(trapz(x, (S_iter - S_intiter).^2  + dt * (dW_dx.^2)));      

                EA_N = [EA_N, err];   
                
                tot_iter_number = tot_iter_number + 1;

                %% Convergence Check
                if err < 1e-10
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
