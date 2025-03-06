%% Biofilm Model - Convergence Order
% This script solves the Biofilm model using different iterative 
% schemes (M scheme, N scheme, M-Adaptive scheme) and visualizes the Convergence order.

clear all;
close all;
clc;

%% Define Parameters
global m_exp;
global u_star;
global phi_ustar;

u_star = 0.367774963378906; 
phi_ustar = 0.0387; 
m_exp = 6;

L = -8;  
R = 8;  
n_values = 100:400:3000;  
T = 1;  dt = 0.1;  ee = 0.5;
Mbb = 0.01;  MBB = 0.01;

%% Solve for Different Methods
[h_M_values, average_orders_M] = M_Order_MethodSolver(L, R, n_values, Mbb, MBB, T, dt, ee, m_exp);

[h_N_values, average_orders_N] = N_Order_MethodSolver(L, R, n_values, T, dt, ee, m_exp);
    
[h_MAdap_values, average_orders_AdapM] = MAdap_Order_MethodSolver(L, R, n_values, T, dt, ee, m_exp);

%% Plotting
hFig = figure;
set(hFig, 'Position', [100, 100, 800, 600]); % Adjust figure position and size

hold on;
grid on;

plot(h_M_values, average_orders_M, 'LineStyle', '-.', 'Color', [0, 0.447, 0.741], ...
    'LineWidth', 5, 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', [0.85, 0.47, 0.32], ...
    'MarkerFaceColor', 'b', 'DisplayName', 'M scheme');

plot(h_MAdap_values, average_orders_AdapM, 'LineStyle', ':', 'Color', [0.466, 0.674, 0.188], ...
    'LineWidth', 5, 'Marker', 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.9, 0.5, 0.77], ...
    'MarkerFaceColor', 'g', 'DisplayName', 'MAdap scheme');

plot(h_N_values, average_orders_N, 'LineStyle', '-.', 'Color', [0.85, 0.325, 0.098], ...
    'LineWidth', 5, 'Marker', '+', 'MarkerSize', 15, 'MarkerEdgeColor', [0.1, 0.25, 0.89], ...
    'MarkerFaceColor', 'r', 'DisplayName', 'N scheme');

xlabel('$1/h$', 'FontSize', 26, 'Interpreter', 'latex');
ylabel('Order of convergence', 'FontSize', 26, 'Interpreter', 'latex');
title('$\tau=0.1$, $T=1$', 'Interpreter', 'latex', 'FontSize', 28);

% Customize grid appearance
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];  % Light gray grid color
ax.GridAlpha = 0.9;  % Grid line transparency
ax.GridLineStyle = '-';  % Solid grid lines
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.Box = 'on';  % Remove the box around the plot

% Set the background color to a light gray
ax.Color = [1, 1, 1];

% Customize legend
legend('Location', 'northeast', 'FontSize', 32, 'Interpreter', 'latex');

% Annotate data points with number of iterations
offset = 0.9; % Vertical offset for the second graph's iteration numbers
 yl = ylim;
 ylim([0 yl(2)+2])
set(gca, 'FontSize', 26);
ax.LineWidth = 2;
% Exporting Figure
exportgraphics(gcf, 'order_BIOFILM.eps', 'BackgroundColor', 'none');
