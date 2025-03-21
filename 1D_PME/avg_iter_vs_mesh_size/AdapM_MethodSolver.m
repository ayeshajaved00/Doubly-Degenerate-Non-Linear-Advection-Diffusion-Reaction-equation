%% PME - average iteration Vs mesh size 
% This script solves the PME using iterative 
% scheme (M adaptive scheme).
function [AdapM_iter_per_step_values, h_MAdap_values] = AdapM_MethodSolver(L, R, n_values, T, dt, ee, m)
    AdapM_iter_per_step_values = zeros(size(n_values));
    h_MAdap_values = zeros(size(n_values));
    Mb = 1;
    MB = 1;

    %% Adaptivity parameters
    adap_num_values = 297;
    adap_start_value = -10;
    adap_end_value = -2;
    adap_step_size = (adap_end_value - adap_start_value) / (adap_num_values - 1);

    for n_index = 1:length(n_values)
        n = n_values(n_index);
        dx = (R - L) / n;
        h = 1 / dx;
        x = L + dx * (0:n);
        
        %% Number of Time Steps
        N = floor(T / dt) + 1;
        t = zeros(1, N);
        
        %% Initial and Boundary Conditions
        Pl = 0; Pr = 0;
        S0 = InitialWW(x, m, 0, 1);
        S_current = S0;
        b_n = S0;
        
        %% Iteration Setup
        iteration_no = 50;
        diverge = 0;
        tot_iter_number = 0;

        %% Precompute Matrices
        tm = (dt / dx^2) .* ones(1, n - 1);
        tp = (dt / dx^2) .* ones(1, n - 1);
        Main_D = [1, (tm + tp), 1];
        Lower_D = [0, -tm, 0];
        Upper_D = [0, -tp, 0];
        
        for j = 1:N
            t(j) = (j - 1) * dt;
            EA = [];
            S_iter = S_current;
            Mb(1) = Mb(end); Mb = [Mb(1)];
            MB(1) = MB(end); MB = [MB(1)];
            
            %% Compute b, B, and their derivatives
            b_iter = bs(S_iter, m);
            B_iter = BBs(S_iter, m);
            dBds_iter = dBBds(S_iter, m);
            dbds_iter = dbds(S_iter, m);
            
            for i = 1:iteration_no
                %% Compute Lb and LB
                Lb = min(max(dbds_iter + Mb(i) * dt, 2 * Mb(i) * dt), 1);
                LB = min(max(dBds_iter + MB(i) * dt, 2 * MB(i) * dt), 1);
                
                AA = Main_D + (Lb ./ LB) .* (1 - ee * dt);
                AA(1) = 1; AA(end) = 1;
                
                f = b_n(2:n) - (1 - ee * dt) .* b_iter(2:n) + (1 - ee * dt) .* B_iter(2:n) .* (Lb(2:n) ./ LB(2:n));
                F = [(1 + (LB(1) / Lb(1))) * Pl, f, (1 + (LB(end) / Lb(end))) * Pr];
                
                if i == 1
                    W_intiter = B_iter;
                else
                    W_intiter = W_iter;
                end
                W_iter = tridiagQQ(AA, Lower_D, Upper_D, F)';
                
                S_intiter = S_iter;
                S_iter = (1 ./ LB) .* (W_iter - B_iter) + S_intiter;
                U_iter = b_iter + Lb .* (S_iter - S_intiter);
                
                %% Update b, B, and their derivatives
                b_iter = bs(S_iter, m);
                B_iter = BBs(S_iter, m);
                dBds_iter = dBBds(S_iter, m);
                dbds_iter = dbds(S_iter, m);
                
                %% Error Calculation
                dW_dx = gradient(W_iter - W_intiter, dx);
                err = sqrt(trapz(x, (Lb .* LB .* (S_iter - S_intiter).^2 + dt *dW_dx.^2)));
                EA = [EA, err];
                
                tot_iter_number = tot_iter_number + 1;
                
                %% Convergence Check
                if err < 1e-6
                    break;
                elseif isnan(err) || err > 100
                    diverge = 1;
                    break;
                end
                
                %% Adaptivity
                if rem(i, 3) == 2
                    U = U_iter - b_iter;
                    O = W_iter - B_iter;
                    
                    for j1 = 0:adap_num_values
                        value = 10^(adap_start_value + j1 * adap_step_size);
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
            end
            
            %% Update Current Solution
            S_current = S_iter;
            b_n = bs(S_current, m);
            
            if diverge == 1
                iter_per_timestep = inf;
                break;
            end
        end
        
        if diverge ~= 1
            iter_per_timestep = tot_iter_number / N;
        end
        
        AdapM_iter_per_step_values(n_index) = iter_per_timestep;
        h_MAdap_values(n_index) = h;
    end
end
