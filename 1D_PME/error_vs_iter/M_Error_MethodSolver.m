function [M_iter_per_step_values, h_M_values, EA_M] = M_Error_MethodSolver(L, R, n_values, T, dt, ee, m)
    
    %% Initialize arrays to store results
    M_iter_per_step_values = zeros(size(n_values));
    h_M_values = zeros(size(n_values));
    err_M_values = zeros(size(n_values));
    iteration_M_iter = zeros(size(err_M_values));
    Mb = 0.001; 
    MB = 0.001;
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
        b_n = S0;
        S_current = S0;

        %% Basic iteration information
        t_start = 1; 
        t_end = N;

        EA_M = []; 
        errors = [];
        iteration_no = 50;
        diverge = 0;

        %% Start of time loop
        tot_iter_number = 0;

        for j = t_start:t_end
            t(j) = (j - 1) * dt;

            if rem(j, 1) == 0
                EA_M = []; 
                errors = [];
            end

            S_iter = S_current;

            %% Precomputing the arrays
            tm = ((dt / dx^2)) .* ones(1, n - 1);
            tp = ((dt / dx^2)) .* ones(1, n - 1);

            Main_D = [0, (tm + tp), 0];
            Lower_D = [0, -tm, 0];
            Upper_D = [0, -tp, 0];

            for i = 1:iteration_no
                %% Update b and B
                b_iter = bs(S_iter, m); % This is the value of b at S^{i-1}_n
                B_iter = BBs(S_iter, m); % This is the value of B at S^{i-1}_n
                dBds_iter = dBBds(S_iter, m); % This is the value of B' at S^{i-1}_n
                dbds_iter = dbds(S_iter, m); % This is the value of b' at S^{i-1}_n

                %% Definition $L_b$ and $L_B$
                Lb = min(max(dbds_iter + Mb * dt, 2 * Mb * dt), 1);
                LB = min(max(dBds_iter + MB * dt, 2 * MB * dt), 1);

                AA = Main_D + (Lb ./ LB) .* (1 - ee * dt);
                AA(1) = 1; 
                AA(end) = 1;

                f = b_n(2:n) - (1 - ee * dt) .* b_iter(2:n) + (1 - ee * dt) .* B_iter(2:n) .* ((Lb(2:n)) ./ (LB(2:n)));
                F = [(1 + (LB(1) / Lb(1))) .* Pl, f, (1 + (LB(end) / Lb(end))) .* Pr];

                if i == 1
                    W_intiter = B_iter;
                else
                    W_intiter = W_iter;
                end

                W_iter = tridiagQQ(AA, Lower_D, Upper_D, F)';
                S_intiter = S_iter;
                S_iter = (1 ./ LB) .* (W_iter - B_iter) + S_intiter;

                %% Error Calculation
                dW_dx = gradient(W_iter - W_intiter, dx);
                
                err = sqrt(trapz(x, (S_iter - S_intiter).^2 + dt * (dW_dx.^2)));

                EA_M = [EA_M, err];
                
                tot_iter_number = tot_iter_number + 1;

                %% Convergence Check
                if err < 1e-6
                    break;
                elseif isnan(err) || err > 100
                    diverge = 1;
                    break;
                end
            end

            %% Variable Update
            S_current = S_iter;
            b_n = bs(S_current, m);

            if diverge == 1
                iter_per_timestep = inf;
                break;
            end
        end

        if diverge ~= 1
            iter_per_timestep = tot_iter_number / t_end;
        end

        M_iter_per_step_values(n_index) = iter_per_timestep;
        h_M_values(n_index) = h;
        iteration_M_iter(n_index) = i;
    end
end