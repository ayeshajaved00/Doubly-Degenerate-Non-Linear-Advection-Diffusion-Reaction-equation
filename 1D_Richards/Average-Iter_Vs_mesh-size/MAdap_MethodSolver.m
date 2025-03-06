%% Richards Equation - Average iteration Vs mesh size (h)
% This script solves the Richards equation using iterative 
% scheme (M-Adaptive scheme) and visualizes the average iteration Vs mesh size (h)

function [MAdap_iter_per_step_values, h_MAdap_values] = MAdap_MethodSolver(L, R, n_values, b_values, B_values, T, dt, ee, m)
   
    global Ks
    global lambda
    
    %% Initialize arrays to store data

    MAdap_iter_per_step_values = zeros(size(n_values));
    h_MAdap_values = zeros(size(n_values));

    %% Adaptivity parameters
    Mb = 1; MB = 1;
    adap_num_values = 16; % Number of values in the range
    adap_start_value = -10; % Starting value
    adap_end_value = -2; % Ending value
    adap_step_size = (adap_end_value - adap_start_value) / (adap_num_values -1);

    for n_index = 1:length(n_values)
        n = n_values(n_index);
        dx = (R-L)/n;
        h = 1 / dx;
        x = zeros(1, n+1) + 1;
        for k = 1:n+1
            x(k) = L + dx * (k - 1);
        end

        %% The Number of Time Steps
        N = floor(T/dt) + 1;
        t = zeros(1, N);
        
        %% Boundary Condition
        Pl = 0;
        Pr = 0;
        
        %% Initial condition
        S0 = InitialWW(x, m, 0, 1);
        b_n = S0;
        S_current = S0;
        
        %% Basic iteration information
        t_start = 1; t_end = N;
        diverge = 0;
        iteration_no = 50;
        tot_iter_number = 0;
        
        %% Start of time loop
        for j = t_start:t_end 
            t(j) = (j - 1) * dt;
            EA_MAdap = []; 
            S_iter = S_current;
            
            %% Precomputing the arrays
            tm = ((dt/dx^2)) .* ones(1, n-1);
            tp = ((dt/dx^2)) .* ones(1, n-1);
            Main_D = [1 (tm + tp) 1];
            Lower_D = [0 -tm 0];
            Upper_D = [0 -tp 0];
            Mb(1) = Mb(end); Mb = [Mb(1)];
            MB(1) = MB(end); MB = [MB(1)];
            
            %% Initial variable update
            S_iter = max(S_iter, 0);
            index = idx(S_iter, b_values);
            b_iter = bs(S_iter, index, b_values);
            B_iter = BB(S_iter, index, B_values);
            dbds_iter = dbds(S_iter, index, b_values);
            dBds_iter = dBBds(S_iter, index, B_values);
            F_iter = Ks(b_iter, lambda);
            
            for i = 1:iteration_no
                %% Defining LB and Lb
                Lb = min(max(dbds_iter + Mb(i) * dt, 2 * Mb(i) * dt), 1);
                LB = min(max(dBds_iter + MB(i) * dt, 2 * MB(i) * dt), 1);
                
                AA = Main_D + (Lb ./ LB) .* (1 - ee * dt);
                AA(1) = 1; AA(end) = 1;
                f = b_n(2:n) - (1 - ee * dt) .* b_iter(2:n) + (1 - ee * dt) .* B_iter(2:n) .* ((Lb(2:n)) ./ (LB(2:n))) + (dt/dx) .* (F_iter(1:n-1) - F_iter(2:n));
                F = [Pl f Pr];
                
                if i == 1
                   W_intiter = B_iter;
                else
                   W_intiter = W_iter;
                end
                
                W_iter = tridiagQQ(AA, Lower_D, Upper_D, F)';
                
                S_intiter = S_iter;
                F_intiter = F_iter;
                S_iter = (1 ./ LB) .* (W_iter - B_iter) + S_intiter;
                U_iter = b_iter + Lb .* (S_iter - S_intiter);
                S_iter = max(S_iter, 0);
                index = idx(S_iter, b_values);
                b_iter = bs(S_iter, index, b_values);
                B_iter = BB(S_iter, index, B_values);
                dbds_iter = dbds(S_iter, index, b_values);
                dBds_iter = dBBds(S_iter, index, B_values);
                F_iter = Ks(b_iter, lambda);
                
                %% Error check
                dW_dx = gradient(W_iter - W_intiter, dx);
                                
                err = sqrt(trapz(x, Lb .* LB .* (S_iter - S_intiter).^2 + dt * (dW_dx.^2)));

                EA_MAdap = [EA_MAdap, err];
                
                tot_iter_number = tot_iter_number + 1;
                
                %% Convergence Check
                if err < 1e-6
                    break;
                elseif isnan(err) || (err > 100)
                    diverge = 1;
                    break; 
                end
                
                %% Adaptivity
                if rem(i, 3) == 2
                    U = U_iter - b_iter;
                    O = W_iter - B_iter;
                    dFiter = F_iter - F_intiter;
                    for j1 = 0:adap_num_values
                        value = adap_start_value + (j1) * adap_step_size;
                        value = 10.^value;
                        Lb = min(max(dbds_iter + value * dt, 2 * value * dt), 1);
                        LB = min(max(dBds_iter + value * dt, 2 * value * dt), 1);
                        G = eta(LB./Lb, U, O, dFiter, dt, dx);
                        if (G < err)
                            break;
                        end
                    end
                    Mb(i+1) = value;
                    MB(i+1) = value;
                else
                    Mb(i+1) = Mb(i);
                    MB(i+1) = MB(i);
                end
            end 
            
            %% Variable Update
            S_current = S_iter;
            b_n = bs(S_current, index, b_values);
            if diverge == 1
                iter_per_timestep = inf;
                break;
            end
        end
        
        if diverge ~= 1
            iter_per_timestep = tot_iter_number / t_end;
        end
        MAdap_iter_per_step_values(n_index) = iter_per_timestep;
        h_MAdap_values(n_index) = h;
    end
end
