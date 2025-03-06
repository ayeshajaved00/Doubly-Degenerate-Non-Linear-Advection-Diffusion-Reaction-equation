%% Richards Equation - Contraction rate
% This script solves the Richards equation using different iterative 
% schemes (M-scheme, Newton scheme and M-Adaptive scheme) and visualizes the contraction rate

clc
clear all;
close all;

time_steps = [0.00001,0.000025,0.00005,0.0001,0.00025,0.0005,0.001, 0.0025,0.005, 0.01, 0.025,0.05,0.1,0.25,0.5,1];

%% Data for plotting
Arithmetic_mean_N = [7.056214227563492e-06, 1.600766437918747e-05, 2.249972260084866e-05, 6.810710138593115e-05, ...
                     3.389499487058266e-04, 0.001204525540284, 0.002093295531731, 0.009236313712345, 0.020772108377058, ...
                     0.028425926637399, 0.034873221168899, 0.054868500564138, 0.079532903756246, inf, inf, inf];

Arithmetic_mean_M = [7.714046003566344e-05, 1.924480717448551e-04, 3.883226012040598e-04, 6.189377441797398e-04, ...
                     0.001569841601682, 0.006243833975659, 0.011575200856382, 0.031805360789534, 0.058480036832183, ...
                     0.107175674630735, 0.161962846870766, 0.256822539321787, 0.415782680951759, 0.866527065401595, ...
                     0.961676080236254, 1.026172377858111];

Arithmetic_mean_Madap = [7.716253105006345e-05, 1.924915551216618e-04, 3.884336897945427e-04, 6.194290400391142e-04, ...
                         0.001571757204725, 0.006108036727245, 0.011576802156714, 0.031498591119714, 0.057375546785179, ...
                         0.104604493239346, 0.159556011469969, 0.255630197722322, 0.362189911838963, 0.623865361630491, ...
                         0.908579703907304, 1.026172377858111];

% Create figure
hFig = figure;
set(hFig, 'Position', [100, 100, 800, 600]); % Set figure size and position

% Plot the data
plot(log10(time_steps), log10(Arithmetic_mean_M), 'LineStyle', '-', 'Color', [0, 0.447, 0.741], 'LineWidth', 5, ...
     'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', [0.85 0.47 0.32], 'MarkerFaceColor', 'b', 'DisplayName', 'M scheme');
hold on;
plot(log10(time_steps), log10(Arithmetic_mean_N), 'LineStyle', '-', 'Color', [0.85, 0.325, 0.098], 'LineWidth', 5, ...
     'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', [0.1 0.25 0.89], 'MarkerFaceColor', 'r', 'DisplayName', 'N scheme');
plot(log10(time_steps), log10(Arithmetic_mean_Madap), 'LineStyle', '-.', 'Color', [0.466, 0.674, 0.188], 'LineWidth', 5, ...
     'Marker', 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.9 0.5 0.77], 'MarkerFaceColor', 'g', 'DisplayName', 'MAdap scheme');

% Axis labels and title
xlabel('$\log_{10}(\tau)$', 'FontSize', 28, 'Interpreter', 'latex');
ylabel('$\log_{10}(\alpha)$', 'FontSize', 28, 'Interpreter', 'latex');
title('$h = 0.025$, $T = 1$', 'FontSize', 28, 'Interpreter', 'latex');

% Adjust y-axis limits
ylim([-7 1]);

% Set grid
grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];  % Light gray grid color
ax.GridAlpha = 0.9;  % Grid line transparency
ax.GridLineStyle = '-';  % Solid grid lines
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.Box = 'on';  % Box around the plot

% Set font size for the axes
set(gca, 'FontSize', 24);

% Customize legend
legend('Location', 'southeast', 'FontSize', 26, 'Interpreter', 'latex');

% Set background color and line width
ax.Color = [1, 1, 1];
ax.LineWidth = 1.5;

% Export the figure
exportgraphics(gcf, 'RichardsCont.eps', 'BackgroundColor', 'none');
