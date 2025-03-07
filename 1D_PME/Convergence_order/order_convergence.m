%% PME - Convergence Order
% This script solves the PME using different iterative 
% schemes (M scheme, N scheme, M-Adaptive scheme) and visualizes the convergence order.

clear all;
close all;
clc;

%% Define parameters
L = -8; R = 8;
n_values = 100:3000:12000;
T = 1; dt = 0.1;
ee = 0; m = 6;

%% Solve the problem for different 'n' values
[h_M_values, average_orders_M] = M_Order_MethodSolver(L, R, n_values, T, dt, ee, m);
[h_N_values, average_orders_N] = Newton_Order_MethodSolver(L, R, n_values, T, dt, ee, m);
[h_AdapM_values, average_orders_AdapM] = AdapM_Order_MethodSolver(L, R, n_values, T, dt, ee, m);

%% Plot results
hFig = figure;
set(hFig, 'Position', [100, 100, 800, 600]);

% Plot M scheme
plot(h_M_values, average_orders_M, 'LineStyle', '--', 'Color', [0, 0.447, 0.741], ...
    'LineWidth', 6, 'Marker', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.85, 0.47, 0.32], ...
    'MarkerFaceColor', 'b', 'DisplayName', 'M scheme');
hold on;

% Plot N scheme
plot(h_M_values, average_orders_N, 'LineStyle', '-', 'Color', [0.85, 0.325, 0.098], ...
    'LineWidth', 6, 'Marker', '+', 'MarkerSize', 10, 'MarkerEdgeColor', [0.1, 0.25, 0.89], ...
    'MarkerFaceColor', 'r', 'DisplayName', 'N scheme');

% Plot M-Adap scheme
plot(h_M_values, average_orders_AdapM, 'LineStyle', '--', 'Color', [0.466, 0.674, 0.188], ...
    'LineWidth', 6, 'Marker', 'd', 'MarkerSize', 10, 'MarkerEdgeColor', [0.9, 0.5, 0.77], ...
    'MarkerFaceColor', 'g', 'DisplayName', 'MAdap scheme');

% Label axes
xlabel('$1/h$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('Order of convergence', 'FontSize', 26, 'Interpreter', 'latex');

grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];
ax.GridAlpha = 0.9;
ax.GridLineStyle = '-';
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.Box = 'on';
ax.Color = [1, 1, 1];

% Customize legend
legend('Location', 'northeast', 'FontSize', 35, 'Interpreter', 'latex');

% Set title
title('$\tau = 0.1$, $T = 1$', 'FontSize', 25, 'Interpreter', 'latex');

% Adjust y-axis limits
yl = ylim;
ylim([0, yl(2) + 2]);

% Set font size and line width
set(gca, 'FontSize', 30);
ax.LineWidth = 1.5;

% Export figure
exportgraphics(gcf, 'Porous_order.eps', 'BackgroundColor', 'none');
