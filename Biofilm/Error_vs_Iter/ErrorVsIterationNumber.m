%% Biofilm Model - Error Vs iteration no
% This script shows the plot of Error Vs iter numbers of  the Biofilm model using different iterative 
% schemes (M scheme, N scheme, M-Adaptive scheme).

clear all;
close all;
clc;
%% Define your parameters
global m_exp;
global u_star;
global phi_ustar;

u_star = 0.367774963378906; 
phi_ustar = 0.0387; 
L = -8;  
R = 8;
m_exp=6;

n_values = 1000;

% Define unique colors for each combination of 'l' and scheme
colors_M = [0 0 1; 0 0 0.5]; % Blue shades for M scheme
colors_N = [1 0 0; 0.5 0 0]; % Red shades for N scheme
colors_MAdap = [0 1 0; 0 0.5 0]; % Green shades for MAdap schemeT = 0.1; dt = 0.1;
ee = 0.5;
dt =0.1; T=0.1;

%% Call the function to solve the problem for the array of 'n' values
figure(2);
% Set figure size and position
set(gcf, 'Position', [100, 100, 800, 600]);
hold on;

% Customize legend
legend_info = cell(length(n_values)*3, 1);
legend_index = 1;

for idx = 1:length(n_values)
    n = n_values(idx);
 switch idx
        case 1
            L_value = 10;
 end   


    [M_iter_per_step_values, h_M_values, iteration_M_iter, EA_M,SM_iter] =  M_Error_MethodSolver(L, R, n_values(idx), T, dt, ee, m_exp);
    [N_iter_per_step_values, h_N_values, iteration_N_iter, EA_N,SN_iter] =  N_Error_MethodSolver(L, R, n_values(idx),  T, dt, ee, m_exp);
    [MAdap_iter_per_step_values, h_MAdap_values, iteration_MAdap_iter, EA_MAdap,Mb] = M_Error_Adap_MethodSolver(L, R, n_values(idx), T, dt, ee, m_exp);


    iteration_N_number = 1: length(log10(EA_N));
    iteration_M_number = 1: length(log10(EA_M));
    iteration_MAdap_number = 1: length(log10(EA_MAdap));

     % Plot M scheme
    semilogy(iteration_M_number, log10(EA_M), 'LineStyle', '-', 'Color', colors_M(idx, :), ...
        'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 8);
    legend_info{legend_index} = ['M l=', num2str(L_value)];
    legend_index = legend_index + 1;

    % Plot N scheme
    semilogy(iteration_N_number, log10(EA_N), 'LineStyle', '--', 'Color', colors_N(idx, :), ...
        'LineWidth', 3, 'Marker', 's', 'MarkerSize', 8);
    legend_info{legend_index} = ['N l=', num2str(L_value)];
    legend_index = legend_index + 1;

    % Plot MAdap scheme
    semilogy(iteration_MAdap_number, log10(EA_MAdap), 'LineStyle', ':', 'Color', colors_MAdap(idx, :), ...
        'LineWidth', 3, 'Marker', 'd', 'MarkerSize', 8);
    legend_info{legend_index} = ['MAdap l=', num2str(L_value)];
    legend_index = legend_index + 1;

end

%Plot M_values (corresponding to MAdap l=10)
semilogy(iteration_MAdap_number, log10(Mb), 'LineStyle', '--', 'Color', [0 0.5 0], ...
     'LineWidth', 3, 'Marker', 'd', 'MarkerSize', 8);

legend_info{legend_index} = 'Mval (MAdap l=10) ';
legend_index = legend_index + 1;

grid on;
xlabel('number of iterations', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$\log_{10} (\mathcal{E}_{\mathrm{fix},n})$', 'FontSize', 20, 'Interpreter', 'latex');
legend(legend_info, 'Location', 'northeast', 'FontSize', 24, 'Interpreter', 'latex');
 ylim([-13 5])
 xlim([0 50])
% Customize grid appearance
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];
ax.GridAlpha = 0.5;
ax.Box = 'on';  % Remove the box around the plot
set(gca, 'FontSize', 28, 'LineWidth', 2);
exportgraphics(gcf, 'Bio_ErrVsiter_01.eps', 'BackgroundColor', 'none');
