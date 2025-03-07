%% Biofilm Model - Average iteration Vs Error
% This script solves the Biofilm model using iterative 
% scheme (M adaptive scheme).

function [MAdap_iter_per_step_values, h_MAdap_values] = MAdap_MethodSolver(L, R, n_values, T, dt, ee, m)

    %% Initialize an array to store iter_per_step for each n
    MAdap_iter_per_step_values = zeros(size(n_values));
    h_MAdap_values = zeros(size(n_values));
    
    Mb = 1; 
    MB = 1;

    %% Adaptivity parameters
    adap_num_values = 25;
    adap_start_value = -10;
    adap_end_value = -2;
    adap_step_size = (adap_end_value - adap_start_value) / (adap_num_values - 1);

    for n_index = 1:length(n_values)
        n = n_values(n_index);
        dx = (R - L) / n;
        h = 1 / dx;
        x = linspace(L, R, n + 1);

        %% The Number of Time Steps
        N = floor(T / dt) + 1;
        t = zeros(1, N);

        %% Boundary Condition
        Pl = 0;
        Pr = 0;

        %% Initial Condition
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
            EA_MAdap = [];
            S_iter = S_current;

            %% Precomputing the Arrays
            tm = (dt / dx^2) .* ones(1, n - 1);
            tp = (dt / dx^2) .* ones(1, n - 1);

            Main_D = [1 (tm + tp) 1];
            Lower_D = [0 -tm 0];
            Upper_D = [0 -tp 0];

            Mb(1) = Mb(end);
            Mb = [Mb(1)];
            MB(1) = MB(end);
            MB = [MB(1)];

            %% Basic Variable Update
            b_iter = bs(S_iter);
            B_iter = BBs(S_iter);
            dbds_iter = dbds(S_iter);
            dBds_iter = dBBds(S_iter);

            for i = 1:iteration_no
                %% Defining LB and Lb
                Lb = min(max(dbds_iter + Mb(i) * dt, 2 * Mb(i) * dt), 1);
                LB = min(max(dBds_iter + MB(i) * dt, 2 * MB(i) * dt), 1);

                AA = Main_D + (Lb ./ LB) .* (1 - ee .* (1 - b_n) .* dt);

                f = b_n(2:n) - (1 - ee .* (1 - b_n(2:n)) .* dt) .* b_iter(2:n) + ...
                    (1 - ee .* (1 - b_n(2:n)) .* dt) .* B_iter(2:n) .* ((Lb(2:n)) ./ (LB(2:n)));
                F = [Pl f Pr];
                if i == 1
                   W_intiter = B_iter;
                else
                   W_intiter = W_iter;
                end
                W_iter = tridiagQQ(AA, Lower_D, Upper_D, F)';
                
                S_intiter = S_iter;
                S_iter = (1 ./ LB) .* (W_iter - B_iter) + S_intiter;
                
                U_iter = b_iter + Lb .* (S_iter - S_intiter);
                
                b_iter = bs(S_iter);
                B_iter = BBs(S_iter);
                dbds_iter = dbds(S_iter);
                dBds_iter = dBBds(S_iter);
                
                %% Error Calculation
                dW_dx = gradient(W_iter - W_intiter, dx);
                                
                err = sqrt(trapz(x,  Lb .* LB .*(S_iter - S_intiter).^2 + dt * (dW_dx.^2)));

                EA_MAdap = [EA_MAdap, err];
                
                tot_iter_number = tot_iter_number + 1;
                
                %% Convergence Check
                if err < 1e-6
                    break;
                elseif isnan(err) || err > 100
                    diverge = 1;
                    break;
                end
                
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
                    disp(strcat('M value chosen: ', num2str(Mb(i + 1))));
                else
                    Mb(i + 1) = Mb(i);
                    MB(i + 1) = MB(i);
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
        MAdap_iter_per_step_values(n_index) = iter_per_timestep;
        h_MAdap_values(n_index) = h;
    end
end
