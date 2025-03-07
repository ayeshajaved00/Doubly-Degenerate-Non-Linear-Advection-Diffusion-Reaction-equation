function [N_iter_per_step_values, h_Nx_values] = DDN_MethodSolver(Lx, Rx, Ly, Ry, nx_values, ny_values, T, dt, ee, m)

    N_iter_per_step_values = zeros(1, length(nx_values));
    h_Nx_values = zeros(size(nx_values));

    M = 0; % Initialization of M, which is a constant in your function

    for nx_index = 1:length(nx_values)
        
        nx = nx_values(nx_index);
        ny = ny_values(nx_index); % I assume ny_values should be same length as nx_values

        dx = (Rx - Lx) / nx;
        dy = (Ry - Ly) / ny;

        hx = 1 / dx;

        % Create x and y grids
        x = Lx + dx * (0:nx);  % Using vectorized assignment for x
        y = Ly + dy * (0:ny);  % Using vectorized assignment for y

        % The number of time steps
        N = floor(T / dt) + 1;
        t = zeros(1, N);
        [X, Y] = meshgrid(x, y);

        % Initial Condition
        S0 = DDBarenblatt(X, Y, m, 0, 2);
        S_current = S0;
        b_n = S0;

        % Basic iteration information
        iteration_no = 50;
        diverge = 0;

        % Start of time loop
        tot_iter_number = 0;

        for i = 1:N
            t(i) = (i - 1) * dt;
            S_iter = S_current;
            EA = [];  % Reset EA for each time step

            for V = 1:iteration_no
                B_iter = DDBBs(S_iter);
                dBds_iter = DDdBBds(S_iter);
                b_iter = DDbs(S_iter);
                dbds_iter = DDdbds(S_iter);

                % Definition of Lb and LB
                Lb = min(max(dbds_iter + M * dt, 2 * M * dt), 1);
                LB = min(max(dBds_iter + M * dt, 2 * M * dt), 1);

                % Matrix calculation
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
                        Lmat(ind, ind) = 4 * LB(j, k) * (dt / dx^2) + Lb(j, k) * (1 - ee * dt);
                        Lmat(ind, ind - 1) = -LB(j - 1, k) * (dt / dx^2);
                        Lmat(ind, ind + 1) = -LB(j + 1, k) * (dt / dx^2);
                        Lmat(ind, ind - nx - 1) = -LB(j, k - 1) * (dt / dx^2);
                        Lmat(ind, ind + nx + 1) = -LB(j, k + 1) * (dt / dx^2);

                        f(ind) = b_n(j, k) + LB(j, k) * (2 * dt / dx^2) * S_iter(j, k) - 2 * (dt / dx^2) * B_iter(j, k) - ...
                                 (dt / dx^2) * LB(j - 1, k) * S_iter(j - 1, k) + (dt / dx^2) * B_iter(j - 1, k) - ...
                                 (dt / dx^2) * LB(j + 1, k) * S_iter(j + 1, k) + (dt / dx^2) * B_iter(j + 1, k) - ...
                                 (dt / dy^2) * LB(j, k - 1) * S_iter(j, k - 1) + (dt / dy^2) * B_iter(j, k - 1) - ...
                                 (dt / dy^2) * LB(j, k + 1) * S_iter(j, k + 1) + (dt / dy^2) * B_iter(j, k + 1) + ...
                                 2 * LB(j, k) * (dt / dy^2) * S_iter(j, k) - 2 * (dt / dy^2) * B_iter(j, k) + ...
                                 (1 - ee * dt) * (Lb(j, k) * S_iter(j, k) - b_iter(j, k));
                    end
                end

                A = Lmat;
                S_intiter = S_iter;
                S_iter = A \ f';

                S_iter = reshape(S_iter, size(X));  % Reshape the solution

                if i == 1
                    W_intiter = B_iter;
                else
                    W_intiter = W_iter;
                end

                W_iter = LB .* (S_iter - S_intiter) + B_iter;

                % Error calculation
                [dxW, dyW] = gradient(W_iter - W_intiter, dx, dy);
                err = sqrt(sum((Lb(:) .* LB(:)) .* (S_iter(:) - S_intiter(:)).^2 * dx * dy + dt * dx * dy * (dxW(:).^2 + dyW(:).^2)));

                EA = [EA, err];  % Store error
                tot_iter_number = tot_iter_number + 1;

                % Convergence check
                if err < 1e-6
                    break;
                elseif isnan(err) || err > 10
                    diverge = 1;
                    break;
                end
            end

            % If the simulation diverged, exit the loop
            if diverge == 1
                iter_per_timestep = inf;
                break;
            end

            % Update variables
            S_current = S_iter;
            b_n = DDbs(S_current);
        end

        % Store iteration number and hx values
        if diverge ~= 1
            iter_per_timestep = tot_iter_number / N;
        end

        N_iter_per_step_values(nx_index) = iter_per_timestep;
        h_Nx_values(nx_index) = hx;
    end
end
