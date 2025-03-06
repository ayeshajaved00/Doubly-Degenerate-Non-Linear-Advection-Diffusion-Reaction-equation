%% Clear workspace and initialize
clear all;
close all;
clc;

% Load Richards equation-related functions
Richards;

%% Define parameters
L = -10; R = 10;        % Domain boundaries

n_values = 2000;        % Number of discretization points

% Define color schemes for different methods
colors_M = [0 0 1; 0 0 0.5];       % Blue shades for M scheme
colors_N = [1 0 0; 0.5 0 0];       % Red shades for N scheme
colors_MAdap = [0 1 0; 0 0.5 0];   % Green shades for MAdap scheme

% Time-related parameters
T = 0.1; dt = 0.1;

% Additional parameters
ee = 0.1; m = 6; Mbb = 0.01; MBB = 0.01;

%% Initialize figure for plotting
figure(2);
set(gcf, 'Position', [100, 100, 800, 600]);  % Set figure size
hold on;

% Initialize legend storage
legend_info = cell(length(n_values) * 3, 1);
legend_index = 1;

%% Iterate over different 'n' values
for idx = 1:length(n_values)
    n = n_values(idx);
    
    % Define L_value based on index
    switch idx
        case 1
            L_value = 10;
    end
    
    % Solve the problem for different schemes
    [M_iter_per_step_values, h_M_values, iteration_M_iter, EA_M] = ...
        M_Error_MethodSolver(L, R, n, b_values, B_values, Mbb, MBB, T, dt, ee, m);
    
    [N_iter_per_step_values, h_N_values, iteration_N_iter, EA_N] = ...
        N_Error_MethodSolver(L, R, n, b_values, B_values, T, dt, ee, m);
    
    [MAdap_iter_per_step_values, h_MAdap_values, iteration_MAdap_iter, EA_MAdap, Mb] = ...
        MAdap_Error_MethodSolver(L, R, n, b_values, B_values, T, dt, ee, m);
    
    % Define iteration numbers
    iteration_N_number = 1:length(log10(EA_N));
    iteration_M_number = 1:length(log10(EA_M));
    iteration_MAdap_number = 1:length(log10(EA_MAdap));

    % Plot M scheme
    semilogy(iteration_M_number, log10(EA_M), 'LineStyle', '-', ...
        'Color', colors_M(idx, :), 'LineWidth', 4, 'Marker', 'o', 'MarkerSize', 8);
    legend_info{legend_index} = ['M, l=', num2str(L_value)];
    legend_index = legend_index + 1;

    % Plot N scheme
    semilogy(iteration_N_number, log10(EA_N), 'LineStyle', '--', ...
        'Color', colors_N(idx, :), 'LineWidth', 4, 'Marker', 's', 'MarkerSize', 8);
    legend_info{legend_index} = ['N, l=', num2str(L_value)];
    legend_index = legend_index + 1;

    % Plot MAdap scheme
    semilogy(iteration_MAdap_number, log10(EA_MAdap), 'LineStyle', ':', ...
        'Color', colors_MAdap(idx, :), 'LineWidth', 4, 'Marker', 'd', 'MarkerSize', 8);
    legend_info{legend_index} = ['MAdap, l=', num2str(L_value)];
    legend_index = legend_index + 1;

end

% Plot additional MAdap value
semilogy(1:length(iteration_MAdap_number), log10(Mb), 'LineStyle', '-.', ...
    'Color', colors_MAdap(2, :), 'LineWidth', 3, 'Marker', 'd', 'MarkerSize', 8);
legend_info{legend_index} = 'Mval (MAdap, l=10)';
legend_index = legend_index + 1;

%% Customize plot appearance
grid on;
xlabel('Number of iterations', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$\log_{10} (\mathcal{E}_{\mathrm{fix},n})$', 'FontSize', 20, 'Interpreter', 'latex');
legend(legend_info, 'Location', 'northeast', 'FontSize', 24, 'Interpreter', 'latex');

% Set axis limits
ylim([-12 2]);
xlim([0 50]);

% Customize grid appearance
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];
ax.GridAlpha = 0.5;
ax.Box = 'on';

% Set axis properties
set(gca, 'FontSize', 28, 'LineWidth', 2);

% Export figure
exportgraphics(gcf, 'RichardsErrVsIter_01.eps', 'BackgroundColor', 'none');
