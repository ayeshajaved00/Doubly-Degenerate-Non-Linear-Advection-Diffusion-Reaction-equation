function [N_iter_per_step_values, h_Nx_values, S_iter] = NR_MethodSolver(Lx, Rx, Ly, Ry, nx_values, ny_values, b_values, B_values, T, dt, ee, m)

    global Ks
    global lambda

    N_iter_per_step_values = zeros(1, length(nx_values));
    h_Nx_values = zeros(size(nx_values));
    M = 0;

    for nx_index = 1:length(nx_values)
        nx = nx_values(nx_index);
        ny = ny_values(nx_index);

        dx = (Rx - Lx) / nx;
        dy = (Ry - Ly) / ny;
        hx = 1 / dx;

        x = zeros(1, nx + 1);
        for i = 1:nx + 1
            x(i) = Lx + dx * (i - 1);
        end

        y = zeros(1, ny + 1);
        for i = 1:ny + 1
            y(i) = Ly + dy * (i - 1);
        end
        %% The Number of Time Steps
        N = floor(T / dt) + 1;
        t = zeros(1, N);
        [X, Y] = meshgrid(x, y);

        %% Initial Condition
        S0 = RBarenblatt(X, Y, m, 0, 2);
        S_current = S0;
        b_n = S0;

        %% Basic iteration information
        t_start = 1;
        t_end = N;
        diverge = 0;
        iteration_no = 50;

        %% Start of time loop
        tot_iter_number = 0;

        for i = t_start:t_end
            t(i) = (i - 1) * dt;
            S_iter = S_current;
            EA = [];

            for V = 1:iteration_no
                S_iter = max(S_iter, 0);

                index = idx(S_iter, b_values);

                b_iter = Rbs(S_iter, index, b_values);  % This is the value of b at S^{i-1}_n
                B_iter = RBB(S_iter, index, B_values);  % This is the value of B at S^{i-1}_n
                dbds_iter = Rdbds(S_iter, index, b_values);  % This is the value of b' at S^{i-1}_n
                dBds_iter = RdBBds(S_iter, index, B_values);  % This is the value of B' at S^{i-1}_n

                F_iter = Ks(b_iter, lambda);  % This is the advection term
                dFds_iter = RdFds(b_iter, dbds_iter);
                Fprime_iter = dFds_iter;

                %% Defining LB and Lb (M-scheme)
                Lb = min(max(dbds_iter + M * dt, 2 * M * dt), 1);
                LB = min(max(dBds_iter + M * dt, 2 * M * dt), 1);

                %% matrix calculation
                Lmat = zeros((nx + 1) * (ny + 1), (nx + 1) * (ny + 1));
                f = zeros(1, (nx + 1) * (ny + 1));

                for ind = 1:(nx + 1) * (ny + 1)
                    k = floor(ind / (nx + 1));
                    j = ind - (nx + 1) * k;
                    k = k + 1;

                    if j == 0 || j == 1 || k == 1 || k == ny + 1
                        f(ind) = 0;
                        Lmat(ind, ind) = 1;
                    else
                        Lmat(ind, ind) = 4 .* LB(j, k) .* (dt / dx^2) + Lb(j, k) .* (1 - ee * dt) + 2 * (dt / dx) .* (Fprime_iter(j, k));
                        Lmat(ind, ind - 1) = -LB(j - 1, k) .* (dt / dx^2) - (dt / dx) .* (Fprime_iter(j - 1, k));
                        Lmat(ind, ind + 1) = -LB(j + 1, k) .* (dt / dx^2);
                        Lmat(ind, ind - nx - 1) = -LB(j, k - 1) .* (dt / dx^2) - (dt / dx) .* (Fprime_iter(j, k - 1));
                        Lmat(ind, ind + nx + 1) = -LB(j, k + 1) .* (dt / dx^2);

                        F = b_n(j, k) + LB(j, k) .* (2 * dt / dx^2) .* S_iter(j, k) - 2 * (dt / dx^2) .* B_iter(j, k) - ...
                            (dt / dx^2) .* LB(j - 1, k) .* S_iter(j - 1, k) + (dt / dx^2) .* B_iter(j - 1, k) - ...
                            (dt / dx^2) .* LB(j + 1, k) .* S_iter(j + 1, k) + (dt / dx^2) .* B_iter(j + 1, k) - ...
                            (dt / dy^2) .* LB(j, k - 1) .* S_iter(j, k - 1) + (dt / dy^2) .* B_iter(j, k - 1) - ...
                            (dt / dy^2) .* LB(j, k + 1) .* S_iter(j, k + 1) + (dt / dy^2) .* B_iter(j, k + 1) + ...
                            2 .* LB(j, k) .* (dt / dy^2) .* S_iter(j, k) - 2 .* (dt / dy^2) .* B_iter(j, k) + ...
                            (1 - ee * dt) .* (Lb(j, k) .* S_iter(j, k) - b_iter(j, k)) + ...
                            (dt / dx) .* (-2 .* (F_iter(j, k) - Fprime_iter(j, k) .* S_iter(j, k)) + ...
                            F_iter(j - 1, k) - Fprime_iter(j - 1, k) .* S_iter(j - 1, k) + ...
                            F_iter(j, k - 1) - Fprime_iter(j, k - 1) .* S_iter(j, k - 1));

                        f(ind) = F;
                    end
                end

                A = Lmat;
                S_intiter = S_iter;
                S_iter = A \ f';

                S_iter = reshape(S_iter, size(X));

                if i == 1
                    W_intiter = B_iter;
                else
                    W_intiter = W_iter;
                end
                W_iter = LB .* (S_iter - S_intiter) + B_iter;

                %% error calculation
                [dxW, dyW] = gradient(W_iter - W_intiter, dx, dy);
                err = sqrt(sum((Lb(:) .* LB(:)) .* (S_iter(:) - S_intiter(:)).^2 * dx * dy + dt * dx * dy * (dxW(:).^2 + dyW(:).^2)));

                EA = [EA, err];
                tot_iter_number = tot_iter_number + 1;

                %% Convergence check
                if err < 1e-6
                    break;
                elseif isnan(err) || err > 10
                    diverge = 1;
                    break;
                end
            end

            %% Variable Update
            S_current = S_iter;
            b_n = Rbs(S_current, index, b_values);

            EA;

            % Check if the simulation diverged
            if diverge == 1
                iter_per_timestep = inf;
                break;
            end
        end

        if diverge ~= 1
            iter_per_timestep = tot_iter_number / t_end;
        end
        N_iter_per_step_values(nx_index) = iter_per_timestep;
        h_Nx_values(nx_index) = hx;

    end
end
