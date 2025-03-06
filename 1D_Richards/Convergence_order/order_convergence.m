%% Richards Equation - Order of Convergence Analysis
% This script solves the Richards equation using different iterative 
% schemes (M scheme, N scheme, M-Adaptive scheme) and visualizes 
% the order of convergence.

% Clear workspace and close figures
clearvars; 
close all;
clc;

% Load necessary parameters and functions
Richards;  % Ensure Richards.m is in the same directory or MATLAB path

%% Define Parameters
L = -10;                % Left boundary
R = 10;                 % Right boundary
n_values = 100:500:4000; % Array of grid sizes

T = 1;                  % Final time
dt = 0.1;               % Time step
ee = 0.1;               % Error tolerance
m = 6;                  % Parameter for the method

Mbb = 0.01;             % Parameter for M scheme
MBB = 0.01;             % Parameter for M scheme

%% Solve the problem for different 'n' values

% M scheme
[h_M_values, average_orders_M] = M_Order_MethodSolver(L, R, n_values, ...
    b_values, B_values, Mbb, MBB, T, dt, ee, m);

% N scheme
[h_N_values, average_orders_N] = N_Order_MethodSolver(L, R, n_values, ...
    b_values, B_values, T, dt, ee, m);

% M-Adaptive scheme
[h_MAdap_values, average_orders_AdapM, orders] = MAdap_Order_MethodSolver(L, R, ...
    n_values, b_values, B_values, T, dt, ee, m);

%% Plotting Results

% Create figure
hFig = figure;
set(hFig, 'Position', [100, 100, 800, 600]); % Set figure size and position

% Plot M scheme
plot(h_M_values, average_orders_M, 'LineStyle', '-', 'Color', [0.1, 0.05, 0.58], ...
    'LineWidth', 5, 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', [0.5, 0.85, 0.91], ...
    'MarkerFaceColor', 'b', 'DisplayName', 'M scheme');
hold on;

% Plot N scheme
plot(h_M_values, average_orders_N, 'LineStyle', '-', 'Color', [0.85, 0.325, 0.098], ...
    'LineWidth', 5, 'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', [0.1, 0.25, 0.89], ...
    'MarkerFaceColor', 'r', 'DisplayName', 'N scheme');

% Plot M-Adaptive scheme
plot(h_M_values, average_orders_AdapM, 'LineStyle', '--', 'Color', [0.466, 0.674, 0.188], ...
    'LineWidth', 5, 'Marker', 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.9, 0.5, 0.77], ...
    'MarkerFaceColor', 'g', 'DisplayName', 'MAdap scheme');

% Labels and title
xlabel('$1/h$', 'FontSize', 28, 'Interpreter', 'latex');
ylabel('Order of convergence', 'FontSize', 28, 'Interpreter', 'latex');
title('$\tau=0.1$, $T=1$', 'Interpreter', 'latex', 'FontSize', 28);

% Grid customization
grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];  % Light gray grid color
ax.GridAlpha = 0.9;              % Grid line transparency
ax.GridLineStyle = '-';          % Solid grid lines
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.Box = 'on';                   % Keep plot box

% Background color
ax.Color = [1, 1, 1];

% Customize legend
legend('Location', 'northeast', 'FontSize', 35, 'Interpreter', 'latex');

% Adjust y-axis limits
yl = ylim;
ylim([0 yl(2) + 2]);

% Set font size and line width
set(gca, 'FontSize', 30);
ax.LineWidth = 1.5;

exportgraphics(gcf, 'orderRichards.eps', 'BackgroundColor', 'none');
