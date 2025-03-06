function [MAdap_iter_per_step_values, h_MAdap_values, iteration_MAdap_iter, EA_MAdap, Mb] = ...
         MAdap_Error_MethodSolver(L, R, n_values, b_values, B_values, T, dt, ee, m)
    
    global Ks lambda
    
    MAdap_iter_per_step_values = zeros(size(n_values));
    h_MAdap_values = zeros(size(n_values));
    iteration_MAdap_iter = zeros(size(n_values));
    Mb = 1;
    MB = 1;
    
    %% Adaptivity parameters
    adap_num_values = 295;
    adap_start_value = -10;
    adap_end_value = -2;
    adap_step_size = (adap_end_value - adap_start_value) / (adap_num_values - 1);
    
    for n_index = 1:length(n_values)
        n = n_values(n_index);
        dx = (R - L) / n;
        h = 1 / dx;
        x = L + dx * (0:n);
        
        %% Time parameters
        N = floor(T / dt) + 1;
        t = zeros(1, N);
        
        %% Boundary conditions
        Pl = 0; Pr = 0;
        
        %% Initial condition
        S0 = InitialWW(x, m, 0, 1);
        S_current = S0;
        b_n = S0;
        
        %% Iteration parameters
        EA_MAdap = [];
        diverge = 0;
        iteration_no = 50;
        tot_iter_number = 0;
        
        %% Time loop
        for j = 1:N
            t(j) = (j - 1) * dt;
            if rem(j, 1) == 0, EA_MAdap = []; end
            S_iter = max(S_current, 0);
            
            %% Precompute arrays
            tm = (dt / dx^2) * ones(1, n-1);
            tp = (dt / dx^2) * ones(1, n-1);
            Main_D = [1, (tm + tp), 1];
            Lower_D = [0, -tm, 0];
            Upper_D = [0, -tp, 0];
            
            Mb(1) = Mb(end); MB(1) = MB(end);
            Mb = [Mb(1)]; MB = [MB(1)];
            
            %% Compute required values
            index = idx(S_iter, b_values);
            b_iter = bs(S_iter, index, b_values);
            B_iter = BB(S_iter, index, B_values);
            dbds_iter = dbds(S_iter, index, b_values);
            dBds_iter = dBBds(S_iter, index, B_values);
            F_iter = Ks(b_iter, lambda);
            
            for i = 1:iteration_no
                Lb = min(max(dbds_iter + Mb(i) * dt, 2 * Mb(i) * dt), 1);
                LB = min(max(dBds_iter + MB(i) * dt, 2 * MB(i) * dt), 1);
                
                AA = Main_D + (Lb ./ LB) .* (1 - ee * dt);
                AA(1) = 1; AA(end) = 1;
                
                f = b_n(2:n) - (1 - ee * dt) .* b_iter(2:n) + ...
                    (1 - ee * dt) .* B_iter(2:n) .* (Lb(2:n) ./ LB(2:n)) + ...
                    (dt / dx) .* (F_iter(1:n-1) - F_iter(2:n));
                F = [Pl, f, Pr];
                
                W_intiter = B_iter;
                if i > 1, W_intiter = W_iter; end
                W_iter = tridiagQQ(AA, Lower_D, Upper_D, F)';
                
                S_intiter = S_iter;
                F_intiter = F_iter;
                S_iter = (1 ./ LB) .* (W_iter - B_iter) + S_intiter;
                U_iter = b_iter + Lb .* (S_iter - S_intiter);
                S_iter = max(S_iter, 0);
                
                %% Update values
                index = idx(S_iter, b_values);
                b_iter = bs(S_iter, index, b_values);
                B_iter = BB(S_iter, index, B_values);
                dbds_iter = dbds(S_iter, index, b_values);
                dBds_iter = dBBds(S_iter, index, B_values);
                F_iter = Ks(b_iter, lambda);
                
                %% Error check
                dW_dx = gradient(W_iter - W_intiter, dx);
                err = sqrt(trapz(x, (S_iter - S_intiter).^2 + dt * (dW_dx.^2)));
                EA_MAdap = [EA_MAdap, err];
                tot_iter_number = tot_iter_number + 1;
                
                %% Convergence check
                if err < 1e-10, break; end
                if isnan(err) || err > 100, diverge = 1; break; end
                
                %% Adaptivity
                if rem(i, 3) == 2
                    U = U_iter - b_iter;
                    O = W_iter - B_iter;
                    dFiter = F_iter - F_intiter;
                    
                    for j1 = 0:adap_num_values
                        value = 10^(adap_start_value + j1 * adap_step_size);
                        Lb = min(max(dbds_iter + value * dt, 2 * value * dt), 1);
                        LB = min(max(dBds_iter + value * dt, 2 * value * dt), 1);
                        G = eta(LB ./ Lb, U, O, dFiter, dt, dx);
                        if G < err, break; end
                    end
                    Mb(i+1) = value; MB(i+1) = value;
                else
                    Mb(i+1) = Mb(i); MB(i+1) = MB(i);
                end
            end
            
            %% Update variables
            S_current = S_iter;
            b_n = bs(S_current, index, b_values);
            if diverge, iter_per_timestep = inf; break; end
        end
        
        %% Store iteration results
        if ~diverge, iter_per_timestep = tot_iter_number / N; end
        MAdap_iter_per_step_values(n_index) = iter_per_timestep;
        h_MAdap_values(n_index) = h;
        iteration_MAdap_iter(n_index) = i;
    end
end