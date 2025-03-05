clear all;
close all;
clc;

Richards;

%% Define your parameters
L = -10;  R = 10;

n_values =100:500:8000; 

T = 1; dt = 0.01; 

ee = 0.1; m = 6;
%% Call the function to solve the problem for the array of 'n' values

[M_iter_per_step_values, h_M_values] = M_MethodSolver(L, R, n_values,b_values,B_values, T, dt,ee,m);

[N_iter_per_step_values, h_N_values,EA_N] = N_MethodSolver(L, R, n_values,b_values,B_values, T, dt,ee,m);

[MAdap_iter_per_step_values, h_MAdap_values,Mb] = MAdap_MethodSolver(L, R, n_values, b_values,B_values,T, dt,ee,m);


% Plot M_iter_per_step_values
figure;
% Set figure size and position
set(gcf, 'Position', [100, 100, 800, 600]);

% Plot M_iter_per_step_values with custom style
plot(h_M_values, M_iter_per_step_values, 'LineStyle', '-', 'Color', [0, 0.447, 0.741], 'LineWidth', 4, 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', [0.85 0.47 0.32], 'MarkerFaceColor', 'b', 'DisplayName', 'M scheme');
hold on;

% Plot N_iter_per_step_values with custom style
plot(h_N_values, N_iter_per_step_values, 'LineStyle', '-', 'Color', [0.85, 0.325, 0.098], 'LineWidth', 4, 'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', [0.1 0.25 0.89], 'MarkerFaceColor', 'r', 'DisplayName', 'N scheme');

% Plot MAdap_iter_per_step_values with custom style
plot(h_MAdap_values, MAdap_iter_per_step_values, 'LineStyle', '-', 'Color', [0.466, 0.674, 0.188], 'LineWidth', 4, 'Marker', 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.9 0.5 0.77], 'MarkerFaceColor', 'g', 'DisplayName', 'MAdap scheme');

% Add labels and title with latex interpretation
xlabel('$1/{h}$', 'FontSize', 36, 'Interpreter', 'latex');
ylabel('Avg iter per step', 'FontSize', 36, 'Interpreter', 'latex');
%title('$m = 6$ \& $C = 0.1$', 'interpreter', 'latex', 'FontSize', 18);

% Display the grid
grid on;

% Customize grid appearance
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];  % Light gray grid color
ax.GridAlpha = 0.5;  % Grid line transparency
ax.GridLineStyle = '-';  % Solid grid lines
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.Box = 'on';  % Remove the box around the plot

% Set the background color to a light gray
ax.Color = [1, 1, 1];

% Customize legend
legend('Location', 'northwest', 'FontSize', 38, 'Interpreter', 'latex');

offset = 0.8; % Vertical offset for the second graph's iteration numbers

yl = ylim;
ylim([0 30])

% Increase the font size of axis labels and tick labels
set(gca, 'FontSize', 33);

% Make the axis and grid lines more prominent
ax.LineWidth = 1.5;

% Save the plot as an image with a transparent background
exportgraphics(gcf, 'richards0.01.eps', 'BackgroundColor', 'none');
