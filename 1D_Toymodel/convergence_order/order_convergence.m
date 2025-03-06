%% DD Toy-model - Convergence Order
% This script solves the DD Toy model using different iterative 
% schemes (M scheme, N scheme, M-Adaptive scheme) and visualizes 
% the Convergence order.

clear all;
close all;
clc;

%% Define Parameters
L = -10;  % Left boundary
R = 10;   % Right boundary

n_values = 300:1000:8000;  % Array of mesh sizes (n values)

T = 1;    dt = 0.1;  ee = 0.5;  m = 6;     

% Parameter for M scheme
Mbb = 0.01;  
MBB = 0.01;  

%% Call the functions to solve the problem for the array of 'n' values

% Solve for M scheme
[h_M_values, average_orders_M] = M_Order_MethodSolver(L, R, n_values, Mbb, MBB, T, dt, ee, m);

% Solve for N scheme
[h_N_values, average_orders_N] = Newton_Order_MethodSolver(L, R, n_values, T, dt, ee, m);

% Solve for M-Adaptive scheme
[h_Adap_values, average_orders_AdapM] = AdapM_Order_MethodSolver(L, R, n_values, T, dt, ee, m);

%% Plotting the results

% Create a figure for plotting
hFig = figure;
set(hFig, 'Position', [100, 100, 800, 600]); % Adjust position and size of the figure

% Plot M scheme convergence order 
plot(h_M_values, average_orders_M, 'LineStyle', '-.', 'Color', [0, 0.447, 0.741], ...
    'LineWidth', 6, 'Marker', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.85 0.47 0.32], ...
    'MarkerFaceColor', 'b', 'DisplayName', 'M scheme');
hold on;

% Plot M-Adaptive scheme convergence order
plot(h_N_values, average_orders_AdapM, 'LineStyle', '-', 'Color', [0.466, 0.674, 0.188], ...
    'LineWidth', 6, 'Marker', 'd', 'MarkerSize', 10, 'MarkerEdgeColor', [0.9 0.5 0.77], ...
    'MarkerFaceColor', 'g', 'DisplayName', 'MAdap scheme');

% Plot N scheme convergence order 
plot(h_N_values, average_orders_N, 'LineStyle', '-.', 'Color', [0.85, 0.325, 0.098], ...
    'LineWidth', 6, 'Marker', 's', 'MarkerSize', 10, 'MarkerEdgeColor', [0.1 0.25 0.89], ...
    'MarkerFaceColor', 'r', 'DisplayName', 'N scheme');

% Add labels and title
xlabel('$1/h$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('Order of convergence', 'FontSize', 26, 'Interpreter', 'latex');
title('$\tau=0.1$, $T=1$', 'Interpreter', 'latex', 'FontSize', 28);

% Enable grid and customize appearance
grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];  % Light gray grid color
ax.GridAlpha = 0.9;  % Grid line transparency
ax.GridLineStyle = '-';  % Solid grid lines
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.Box = 'on';  % Remove the box around the plot

% Set background color
ax.Color = [1, 1, 1];

% Customize legend appearance
legend('Location', 'northeast', 'FontSize', 36, 'Interpreter', 'latex');

% Adjust y-axis limits
yl = ylim;
ylim([0 yl(2) + 2]);

% Increase font size for axis labels and tick labels
set(gca, 'FontSize', 25);

% Make the axis and grid lines more prominent
ax.LineWidth = 1.5;

% Save the plot as an EPS file with a transparent background
exportgraphics(gcf, 'DD_Conv_order.eps', 'BackgroundColor', 'none');
