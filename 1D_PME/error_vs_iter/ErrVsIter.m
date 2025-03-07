clear all;
close all;
clc
%% Define your parameters
L = -8;  
R = 8;
m_exp = 6;
n_values = 10000;

% Define unique colors for each combination of 'l' and scheme
colors_M = [0 0 1; 0 0 0.5]; % Blue shades for M scheme
colors_N = [1 0 0; 0.5 0 0]; % Red shades for N scheme
colors_MAdap = [0 1 0; 0 0.5 0]; % Green shades for MAdap scheme
T = 0.1; dt = 0.1;
ee = 0;

%% Plot settings
figure(2);
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;

% Initialize legend
legend_info = cell(length(n_values) * 3, 1);
legend_index = 1;

for idx = 1:length(n_values)
    n = n_values(idx);
    switch idx
        case 1
            L_value = 100;
    end
    
    [M_iter_per_step_values, h_M_values, EA_M] = M_Error_MethodSolver(L, R, n, T, dt, ee, m_exp);
    [N_iter_per_step_values, h_N_values, EA_N] = N_Error_MethodSolver(L, R, n, T, dt, ee, m_exp);
    [MAdap_iter_per_step_values, h_MAdap_values, EA_MAdap,Mb] = M_Error_Adap_MethodSolver(L, R, n, T, dt, ee, m_exp);

    iteration_N_number = 1:length(log10(EA_N));
    iteration_M_number = 1:length(log10(EA_M));
    iteration_MAdap_number = 1:length(log10(EA_MAdap));

    % Plot M scheme
    semilogy(iteration_M_number, log10(EA_M), 'LineStyle', '-', 'Color', colors_M(idx, :), ...
        'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8);
    legend_info{legend_index} = ['M l=', num2str(L_value)];
    legend_index = legend_index + 1;

    % Plot N scheme
    semilogy(iteration_N_number, log10(EA_N), 'LineStyle', '--', 'Color', colors_N(idx, :), ...
        'LineWidth', 3, 'Marker', 's', 'MarkerSize', 8);
    legend_info{legend_index} = ['N l=', num2str(L_value)];
    legend_index = legend_index + 1;

    % Plot MAdap scheme
    semilogy(iteration_MAdap_number, log10(EA_MAdap), 'LineStyle', ':', 'Color', colors_MAdap(idx, :), ...
        'LineWidth', 3, 'Marker', 'd', 'MarkerSize', 8);
    legend_info{legend_index} = ['MAdap l=', num2str(L_value)];
    legend_index = legend_index + 1;

end

%Plot M_values (corresponding to MAdap)
semilogy(iteration_MAdap_number, log10(Mb), 'LineStyle', '-', 'Color', [0 0.5 0], ...
     'LineWidth', 2, 'Marker', 's', 'MarkerSize', 8);

legend_info{legend_index} = 'Mval (MAdap l=100) ';
legend_index = legend_index + 1;


grid on;
xlabel('number of iterations', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$\log_{10} (\mathcal{E}_{\mathrm{fix},n})$', 'FontSize', 20, 'Interpreter', 'latex');
legend(legend_info, 'Location', 'northeast', 'FontSize', 24, 'Interpreter', 'latex');
% Customize grid appearance
ax = gca;
ax.Box = 'on';  % Remove the box around the plot
ax.GridColor = [0.7, 0.7, 0.7];
ax.GridAlpha = 0.5;
set(gca, 'FontSize', 28, 'LineWidth', 2);
ylim([-10 2])
exportgraphics(gcf, 'PME_ErrvsIter_R01.eps', 'BackgroundColor', 'none');
