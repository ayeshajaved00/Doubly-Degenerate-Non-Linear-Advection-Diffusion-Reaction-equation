function [ M_iter_per_step_values, h_M_x_values] =  DDM_MethodSolver(Lx, Rx, Ly, Ry, nx_values, ny_values, T, dt, ee, m)

    M_iter_per_step_values = zeros(1, length(nx_values));
    h_M_x_values = zeros(size(nx_values));

    Mb = 0.01;
    MB = 0.01;

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

        dtdx = dt / (dx^2);

        %% The Number of Time Steps
        N = floor(T / dt) + 1;  
        t = zeros(1, N);
        [X, Y] = meshgrid(x, y);

        % Initial Condition 
        S0 = DDBarenblatt(X, Y, m, 0, 2);
        S_current = S0;
        b_n = S0;

        %% Basic iteration information
        t_start = 1;
        t_end = N;
        iteration_no = 50; 
        diverge = 0;

        %% Start of time loop  
        tot_iter_number = 0;

        for i = t_start:t_end 
            t(i) = (i - 1) * dt;
            S_iter = S_current;
            EA = [];

            for V = 1:iteration_no
                B_iter = DDBBs(S_iter);
                dBds_iter = DDdBBds(S_iter);
                b_iter = DDbs(S_iter);
                dbds_iter = DDdbds(S_iter); 

                %% Definition of Lb and LB
                % M-scheme
                Lb = min(max(dbds_iter + Mb * dt, 2 * Mb * dt), 1);
                LB = min(max(dBds_iter + MB * dt, 2 * MB * dt), 1);

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
                        Lmat(ind, ind) = (4 * (dtdx)) + (Lb(j, k) ./ LB(j, k)) * (1 - ee * dt); 
                        Lmat(ind, ind - 1) = -(dtdx);
                        Lmat(ind, ind + 1) = -(dtdx);
                        Lmat(ind, ind - nx - 1) = -(dtdx);
                        Lmat(ind, ind + nx + 1) = -(dtdx);

                        f(ind) = b_n(j, k) + (1 - ee * dt) * ((Lb(j, k) ./ LB(j, k)) * B_iter(j, k) - b_iter(j, k));
                    end
                end

                A = Lmat;

                if i == 1
                    W_intiter = B_iter;
                else
                    W_intiter = W_iter;
                end

                W_iter = A \ f'; 
                W_iter = reshape(W_iter, size(X));

                S_intiter = S_iter;
                S_iter = (1 ./ LB) .* (W_iter - B_iter) + S_intiter;
                S_iter = reshape(S_iter, size(X)); 

                %% Error check
                [dxW, dyW] = gradient(W_iter - W_intiter, dx, dy);
                err = sqrt(sum((Lb(:) .* LB(:)) .* (S_iter(:) - S_intiter(:)).^2 * dx * dy + dt * dx * dy * (dxW(:).^2 + dyW(:).^2))); 

                EA = [EA, err];   
                tot_iter_number = tot_iter_number + 1;

                %% Convergence Check
                if err < 1e-6         
                    break;
                elseif isnan(err) || err > 10
                    diverge = 1;
                    break; 
                end
            end 

            %% Variable Update
            S_current = S_iter;
            b_n = DDbs(S_current);

            if diverge == 1
                iter_per_timestep = inf;
                break;
            end
        end

        if diverge ~= 1
            iter_per_timestep = tot_iter_number / t_end;
        end

        M_iter_per_step_values(nx_index) = iter_per_timestep;
        h_M_x_values(nx_index) = hx;
    end
end
