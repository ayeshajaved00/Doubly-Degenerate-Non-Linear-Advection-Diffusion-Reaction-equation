%% DD Toy-model - Contraction rate
% This script solves the DD Toy model using different iterative 
% schemes (M scheme, N scheme, M-Adaptive scheme) and visualizes 
% the contraction rate.

clc;
clear all;
close all;

%% Define Time Steps
time_steps = [0.00001, 0.000025, 0.00005, 0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1];

%% Define Arithmetic Mean Values for Different Schemes
Arithmetic_mean_M = [0.313513292392116, 0.313989799092944, 0.315162834619333, 0.314391479080634, 0.314762716363451, 0.319187259875989, ...
                     0.317377200972413, 0.317921545344219, 0.338169123091373, 0.342060947016903, 0.346380158966785, 0.351745722560142, ...
                     0.365465207579260, 0.877226079911590, 0.783987052306858, 0.288209208988736];

Arithmetic_mean_Adap = [0.002555549378494, 0.001840728434337, 0.003298590380522, 0.005201698021760, 0.006158440333301, 0.016243865155120, ...
                            0.017548436943636, 0.025161830492333, 0.044701503871540, 0.008101132308977, 0.028549166152101, 0.098918285810243, ...
                            0.039167122740015,0.080085014351100,0.056146056400081,0.041386727594425];];

Arithmetic_mean_N = [0.001501420410608,0.003125594727534,0.002612866516749,0.004260708238496,0.007977637164320,0.012742514489587,...
                             0.020325305929531,0.023708367579025,0.044706645395159,0.008100764839659,0.027741864413620,0.110722423350793,...
                             0.039413694319520,0.780437600490464,inf,inf];

%% Create the Plot
hFig = figure;
set(hFig, 'Position', [100, 100, 800, 600]);  % Adjust the position and size as needed

% Plot the M scheme
plot(log10(time_steps), log10(Arithmetic_mean_M), 'LineStyle', '-.', 'Color', [0, 0.447, 0.741], 'LineWidth', 6, 'Marker', 'o', ...
     'MarkerSize', 8, 'MarkerEdgeColor', [0.85 0.47 0.32], 'MarkerFaceColor', 'b', 'DisplayName', 'M scheme');
hold on;

% Plot the N scheme
plot(log10(time_steps), log10(Arithmetic_mean_N), 'LineStyle', '-.', 'Color', [0.85, 0.325, 0.098], 'LineWidth', 6, 'Marker', '+', ...
     'MarkerSize', 8, 'MarkerEdgeColor', [0.1 0.25 0.89], 'MarkerFaceColor', 'r', 'DisplayName', 'N scheme');

% Plot the MAdap scheme
plot(log10(time_steps), log10(Arithmetic_mean_Adap), 'LineStyle', '--', 'Color', [0.466, 0.674, 0.188], 'LineWidth', 6, 'Marker', 'd', ...
     'MarkerSize', 8, 'MarkerEdgeColor', [0.9 0.5 0.77], 'MarkerFaceColor', 'g', 'DisplayName', 'MAdap scheme');

%% Customize Plot Appearance
xlabel('$\log_{10}(\tau)$','FontSize', 20, 'Interpreter', 'latex');
ylabel('$\log_{10}(\alpha)$','FontSize', 20, 'Interpreter', 'latex');

% Enable grid and customize appearance
grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];  % Light gray grid color
ax.GridAlpha = 0.9;  % Grid line transparency
ax.GridLineStyle = '-';  % Solid grid lines
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.Box = 'on';  % Remove the box around the plot

% Set axis limits
ylim([-4 1]);

% Add legend
legend('Location', 'southeast', 'FontSize', 30, 'Interpreter', 'latex');

% Set background color
ax.Color = [1, 1, 1];

% Add title
title('$h = 0.025$, $T = 1$', 'FontSize', 30, 'Interpreter', 'latex');
set(gca, 'FontSize', 25);

% Customize axis line width
ax.LineWidth = 1.5;

%% Save the figure
exportgraphics(gcf, 'DD_Cont_rate.eps', 'BackgroundColor', 'none');
