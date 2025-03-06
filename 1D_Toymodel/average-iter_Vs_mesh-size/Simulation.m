%% DD Toy-model - Average iteration vs Mesh size (h)
% This script solves the DD Toy model using different iterative 
% schemes (M scheme, N scheme, M-Adaptive scheme) and visualizes 
% the Average iteration vs Mesh size (h) behavior.

clear all;
close all;
clc;

%% Define parameters
L = -10;                  % Left boundary
R = 10;                   % Right boundary
n_values = 100:1000:20000; % Range of mesh values
T = 1;   dt = 0.1;  ee = 0.5;   m = 6;                 

%% Solve the problem for different 'n' values using iterative schemes
[M_iter_per_step_values, h_M_values] = M_MethodSolver(L, R, n_values, T, dt, ee, m);
[N_iter_per_step_values, h_N_values] = Newton_MethodSolver(L, R, n_values, T, dt, ee, m);
[AdapM_iter_per_step_values, h_Adap_values] = AdapM_MethodSolver(L, R, n_values, T, dt, ee, m);

%% Plotting the results
figure(1);
set(gcf, 'Position', [100, 100, 800, 600]);

% Plot M scheme
plot(h_M_values, M_iter_per_step_values, ...
    'LineStyle', '-', 'Color', [0, 0.447, 0.741], 'LineWidth', 4, ...
    'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', [0.85, 0.47, 0.32], ...
    'MarkerFaceColor', 'b', 'DisplayName', 'M scheme');
hold on;

% Plot N scheme
plot(h_N_values, N_iter_per_step_values, ...
    'LineStyle', '-', 'Color', [0.85, 0.325, 0.098], 'LineWidth', 4, ...
    'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', [0.1, 0.25, 0.89], ...
    'MarkerFaceColor', 'r', 'DisplayName', 'N scheme');

% Plot MAdap scheme
plot(h_Adap_values, AdapM_iter_per_step_values, ...
    'LineStyle', '--', 'Color', [0.466, 0.674, 0.188], 'LineWidth', 4, ...
    'Marker', 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.9, 0.5, 0.77], ...
    'MarkerFaceColor', 'g', 'DisplayName', 'MAdap scheme');

% Labels and title
xlabel('$1/{h}$', 'FontSize', 36, 'Interpreter', 'latex');
ylabel('Avg iter per step', 'FontSize', 36, 'Interpreter', 'latex');
% title('$m = 6$ \& $C = 0.2$', 'Interpreter', 'latex', 'FontSize', 18);

% Customize grid
grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];  % Light gray grid color
ax.GridAlpha = 0.5;  % Grid transparency
ax.GridLineStyle = '-';  % Solid grid lines
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.Box = 'on';  % Box around the plot

% Set background color to light gray
ax.Color = [1, 1, 1];

% Customize legend
legend('Location', 'northeast', 'FontSize', 35, 'Interpreter', 'latex');

% Set y-axis limits and axis appearance
yl = ylim;
ylim([0, 60]);

% Increase font size for axis labels and tick labels
set(gca, 'FontSize', 30);

% Make axis and grid lines more prominent
ax.LineWidth = 1.5;

% Save the plot as an image with transparent background
exportgraphics(gcf, 'DDouble1D0.1.eps', 'BackgroundColor', 'none');
