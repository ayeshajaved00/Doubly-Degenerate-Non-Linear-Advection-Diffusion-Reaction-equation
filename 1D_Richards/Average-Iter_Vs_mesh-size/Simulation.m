%% Richards Equation - Average iteration Vs mesh size (h)
% This script solves the Richards equation using different iterative 
% schemes (M-scheme, Newton scheme and M-Adaptive scheme) and visualizes
% the average iteration Vs mesh size (h)

clear all;
close all;
clc;

% Call the Richards function for solving the problem
Richards;

%% Define your parameters
L = -10;          % Left boundary of the domain
R = 10;           % Right boundary of the domain

n_values = 100:500:8000;  % Array of 'n' values for the number of intervals

T = 1;             % Total time
dt = 0.1;          % Time step size

ee = 0.1;          % Some parameter (maybe related to the problem setup)
m = 6;             % Another parameter for the method

%% Call the solver functions for different methods
[M_iter_per_step_values, h_M_values] = M_MethodSolver(L, R, n_values, b_values, B_values, T, dt, ee, m);
[N_iter_per_step_values, h_N_values] = N_MethodSolver(L, R, n_values, b_values, B_values, T, dt, ee, m);
[MAdap_iter_per_step_values, h_MAdap_values] = MAdap_MethodSolver(L, R, n_values, b_values, B_values, T, dt, ee, m);

%% Plot results

% Set up the figure window
figure;
set(gcf, 'Position', [100, 100, 800, 600]);

% Plot M_iter_per_step_values vs h_M_values (M scheme)
plot(h_M_values, M_iter_per_step_values, 'LineStyle', '-', 'Color', [0, 0.447, 0.741], 'LineWidth', 4, ...
    'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', [0.85 0.47 0.32], 'MarkerFaceColor', 'b', 'DisplayName', 'M scheme');
hold on;

% Plot MAdap_iter_per_step_values vs h_MAdap_values (MAdap scheme)
plot(h_MAdap_values, MAdap_iter_per_step_values, 'LineStyle', '-', 'Color', [0.466, 0.674, 0.188], 'LineWidth', 4, ...
    'Marker', 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.9 0.5 0.77], 'MarkerFaceColor', 'g', 'DisplayName', 'MAdap scheme');

% Plot N_iter_per_step_values vs h_N_values (N scheme)
plot(h_N_values, N_iter_per_step_values, 'LineStyle', '-', 'Color', [0.85, 0.325, 0.098], 'LineWidth', 4, ...
    'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', [0.1 0.25 0.89], 'MarkerFaceColor', 'r', 'DisplayName', 'N scheme');

% Label the axes with LaTeX formatting
xlabel('$1/{h}$', 'FontSize', 36, 'Interpreter', 'latex');
ylabel('Avg iter per step', 'FontSize', 36, 'Interpreter', 'latex');

% Display the grid
grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];  % Light gray grid color
ax.GridAlpha = 0.5;  % Grid transparency
ax.GridLineStyle = '-';  % Solid grid lines
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.Box = 'on';  % Turn on the box around the plot

% Set the background color to light gray
ax.Color = [1, 1, 1];

% Customize the legend
legend('Location', 'northwest', 'FontSize', 38, 'Interpreter', 'latex');

% Set the limits for the plot
ylim([0 30]);

% Increase font size for axis labels and tick labels
set(gca, 'FontSize', 30);

% Make the axis and grid lines more prominent
ax.LineWidth = 1.5;

% Save the plot as an EPS image with a transparent background
exportgraphics(gcf, 'richards0.1.eps', 'BackgroundColor', 'none');
