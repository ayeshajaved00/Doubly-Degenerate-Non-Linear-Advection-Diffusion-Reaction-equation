%% DD Toy-model - Error Vs Number of Iterations
% This script solves the DD Toy model using different iterative 
% schemes (M scheme, N scheme, M-Adaptive scheme) and visualizes 
% the error vs number of iterations.

clear all;
close all;
clc;

%% Define Parameters
L = -10;          % Left boundary of the domain
R = 10;           % Right boundary of the domain
n_values = 10000; 

T = 0.1;   dt = 0.1;   ee = 0.5;  m = 6;           

% M scheme related parameters
Mbb = 0.01;       % M scheme parameter
MBB = 0.01;       % M scheme parameter

% Color schemes for different methods
colors_M = [0 0 1; 0 0 0.5];      % Blue shades for M scheme
colors_N = [1 0 0; 0.5 0 0];      % Red shades for N scheme
colors_MAdap = [0 1 0; 0 0.5 0];  % Green shades for MAdap scheme

%% Prepare Figure
figure(2);
set(gcf, 'Position', [100, 100, 800, 600]);  % Set figure size
hold on;  % Hold for multiple plots

% Initialize legend information
legend_info = cell(length(n_values) * 3, 1);
legend_index = 1;

%% Loop over n values to calculate and plot errors
for idx = 1:length(n_values)
    n = n_values(idx);
    switch idx
        case 1
            L_value = 50;  
    end

    % Call solver functions for each scheme
    [M_iter_per_step_values, h_M_values, iteration_M_iter, EA_M] = M_Error_MethodSolver(L, R, n, Mbb, MBB, T, dt, ee, m);
    [N_iter_per_step_values, h_N_values, iteration_N_iter, EA_N] = N_Error_MethodSolver(L, R, n, T, dt, ee, m);
    [MAdap_iter_per_step_values, h_MAdap_values, iteration_MAdap_iter, EA_MAdap, Mb] = MAdap_Error_MethodSolver(L, R, n, T, dt, ee, m);

    % Define iteration numbers for plotting (log scale)
    iteration_N_number = 1:length(log10(EA_N));
    iteration_M_number = 1:length(log10(EA_M));
    iteration_MAdap_number = 1:length(log10(EA_MAdap));

    % Gap adjustment for visibility
    gap = 1; % Adjust the gap value as needed

    % Plot M scheme 
    semilogy(iteration_M_number, log10(EA_M), 'LineStyle', '--', 'Color', colors_M(idx, :), ...
        'LineWidth', 4, 'Marker', 'o', 'MarkerSize', 8);
    legend_info{legend_index} = ['M l=', num2str(L_value)];
    legend_index = legend_index + 1;

    % Plot N scheme
    semilogy(iteration_N_number, log10(EA_N), 'LineStyle', '--', 'Color', colors_N(idx, :), ...
        'LineWidth', 4, 'Marker', 's', 'MarkerSize', 8);
    legend_info{legend_index} = ['N l=', num2str(L_value)];
    legend_index = legend_index + 1;

    % Plot MAdap scheme
    semilogy(iteration_MAdap_number, log10(EA_MAdap), 'LineStyle', ':', 'Color', colors_MAdap(idx, :), ...
        'LineWidth', 4, 'Marker', 'd', 'MarkerSize', 8);
    legend_info{legend_index} = ['MAdap l=', num2str(L_value)];
    legend_index = legend_index + 1;

end

semilogy(iteration_MAdap_number, log10(Mb), 'LineStyle', '-', 'Color', [0 0.5 0], ...
    'LineWidth', 3, 'Marker', 'd', 'MarkerSize', 8);
legend_info{legend_index} = 'Mval (MAdap l=50)';
legend_index = legend_index + 1;

%% Grid and Styling
grid on;
xlabel('number of iterations', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$\log_{10} (\mathcal{E}_{\mathrm{fix},n})$', 'FontSize', 20, 'Interpreter', 'latex');

% Display the legend
legend(legend_info, 'Location', 'northeast', 'FontSize', 24, 'Interpreter', 'latex');

% Customize the grid appearance
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];  % Light gray grid color
ax.GridAlpha = 0.5;  % Grid line transparency
ax.GridLineStyle = '-';  % Solid grid lines
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.Box = 'on';  % Box around the plot

% Set background color
ax.Color = [1, 1, 1];
ylim([-13 2]);

% Increase font size and make axis lines more prominent
set(gca, 'FontSize', 28, 'LineWidth', 2);
ax.LineWidth = 2;

%% Save the figure
exportgraphics(gcf, 'DD_ErrvsIter_1e-10_01.eps', 'BackgroundColor', 'none');
