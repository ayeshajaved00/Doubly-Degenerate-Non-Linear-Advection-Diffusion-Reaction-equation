clear all;
close all;
clc;

Richards;

%% Define your parameters
Lx = -10; Rx = 10; Ly = -10; Ry = 10;

T = 0.1; dt = 0.1; ee = 0.1; m = 6;

nx_values = 10:10:30; ny_values = 10:10:30;

%% Call the function to solve the problem for the array of 'n' values
[M_iter_per_step_values, h_M_x_values] = MR_MethodSolver(Lx, Rx, Ly, Ry, nx_values, ny_values, b_values, B_values, T, dt, ee, m);

[N_iter_per_step_values, h_Nx_values, S_iter_N] = NR_MethodSolver(Lx, Rx, Ly, Ry, nx_values, ny_values, b_values, B_values, T, dt, ee, m);

[MAdap_iter_per_step_values, h_MAdap_x_values] = MAdapR_MethodSolver(Lx, Rx, Ly, Ry, nx_values, ny_values, b_values, B_values, T, dt, ee, m);

figure(1);
%% Set figure size and position
set(gcf, 'Position', [100, 100, 800, 600]);

plot(h_M_x_values, M_iter_per_step_values, 'LineStyle', '-', 'Color', [0, 0.447, 0.741], 'LineWidth', 6, 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', [0.85 0.47 0.32], 'MarkerFaceColor', 'b', 'DisplayName', 'M scheme');
hold on;

plot(h_Nx_values, N_iter_per_step_values, 'LineStyle', '--', 'Color', [0.85, 0.325, 0.098], 'LineWidth', 6, 'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', [0.1 0.25 0.89], 'MarkerFaceColor', 'r', 'DisplayName', 'N scheme');

plot(h_MAdap_x_values, MAdap_iter_per_step_values, 'LineStyle', ':', 'Color', [0.466, 0.674, 0.188], 'LineWidth', 6, 'Marker', 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.9 0.5 0.77], 'MarkerFaceColor', 'g', 'DisplayName', 'MAdap scheme');

%% Add labels and title with latex interpretation
xlabel('$1/{h}$', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('Avg iter per step', 'FontSize', 30, 'Interpreter', 'latex');

%% Display the grid
grid on;

%% Customize grid appearance
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];  % Light gray grid color
ax.GridAlpha = 0.5;  % Grid line transparency
ax.GridLineStyle = '-';  % Solid grid lines
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.Box = 'on';  % Remove the box around the plot

%% Set the background color to a light gray
ax.Color = [1, 1, 1];

%% Customize legend
legend('Location', 'northeast', 'FontSize', 36, 'Interpreter', 'latex');

%% Annotate data points with number of iterations
offset = 0.8; % Vertical offset for the second graph's iteration numbers

yl = ylim;
ylim([0 yl(2) + 10])

%% Increase the font size of axis labels and tick labels
set(gca, 'FontSize', 34);

%% Make the axis and grid lines more prominent
ax.LineWidth = 1.5;

%% Store data in excel file 
data = [h_M_x_values.', M_iter_per_step_values.', N_iter_per_step_values.', MAdap_iter_per_step_values.'];
% Excel column headers
headers = {'h_M_x_values', 'MS', 'NS', 'MAdapS'};

% Specify the filename
filename = 'RR01x.xlsx';

% Create a table with the data and headers
T = array2table(data, 'VariableNames', headers);

% Write the table to an Excel file
writetable(T, filename);

%% Save the plot as an image with a transparent background
exportgraphics(gcf, 'Richards2D0.1.eps', 'BackgroundColor', 'none');
