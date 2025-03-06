%% Biofilm Model - Contraction rate
% This script plot the contraction rate of the Biofilm model using different iterative 
% schemes (M scheme, N scheme, M-Adaptive scheme).

clc
clear all;
close all;

time_steps = [0.00001,0.000025,0.00005,0.0001,0.00025,0.0005,0.001, 0.0025,0.005, 0.01, 0.025,0.05,0.1,0.25,0.5,1];

%% h= 0.02 , N=800,  err = sqrt(trapz(x, (S_iter - S_intiter).^2 + dt * (dW_dx.^2)));

 Arithmetic_mean_M = [4.427404744340027e-06,1.520916551476124e-05,1.301758886730194e-04,2.052423838360343e-04,0.001317005199643,0.003000469236908,0.009934534389219,0.042089480033151,0.050580998037892,0.013048953309394,0.053153716691635,0.229543336753767,0.540347646862752,0.809391447751282,0.971641116803704,0.969632860354972];

 Arithmetic_mean_Adap = [4.141017160601342e-06,1.324391842575564e-05,1.018475898923045e-04,4.087492430991211e-04,0.002792220419336,2.820398911653550e-04,6.354507342677017e-04,0.041271369709215,0.049890088533605,0.050260596349608,0.028343150141773,0.057129032768500,0.052918853439674,0.088904938060692,0.056817577069656,0.045998115284744];
 
 Arithmetic_mean_N = [4.122571248079025e-06,1.487909353419424e-05,1.327574379896663e-04,1.888639301216823e-04,0.001759558217307,0.002478340756707,0.010488079795288,0.010207093954378,0.052462767511531,0.050260592481947,0.065231393281297,0.053904826208888,0.074483482523228,0.047945699141989,inf,inf];

hFig = figure;
set(hFig, 'Position', [100, 100, 900, 700]); % Slightly larger figure for clarity

% Define colors
color_M = [0, 0.447, 0.741];
color_N = [0.85, 0.325, 0.098];
color_Adap = [0.466, 0.674, 0.188];

% Plot M scheme
plot(log10(time_steps), log10(Arithmetic_mean_M), ...
    'LineStyle', ':', 'Color', color_M, 'LineWidth', 4, ...
    'Marker', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.85 0.47 0.32], ...
    'MarkerFaceColor', 'b', 'DisplayName', 'M scheme');
hold on;

% Plot N scheme
plot(log10(time_steps), log10(Arithmetic_mean_N), ...
    'LineStyle', '-', 'Color', color_N, 'LineWidth', 4, ...
    'Marker', 's', 'MarkerSize', 10, 'MarkerEdgeColor', [0.1 0.25 0.89], ...
    'MarkerFaceColor', 'r', 'DisplayName', 'N scheme');

% Plot Adaptive M scheme
plot(log10(time_steps), log10(Arithmetic_mean_Adap), ...
    'LineStyle', '-.', 'Color', color_Adap, 'LineWidth', 4, ...
    'Marker', 'd', 'MarkerSize', 10, 'MarkerEdgeColor', [0.9 0.5 0.77], ...
    'MarkerFaceColor', 'g', 'DisplayName', 'MAdap scheme');

% Set axis labels
xlabel('$\log_{10}(\tau)$', 'FontSize', 32, 'Interpreter', 'latex');
ylabel('$\log_{10}(\alpha)$', 'FontSize', 32, 'Interpreter', 'latex');

% Adjust axes limits
ylim([-6 2]); % Set proper limits
xlim([min(log10(time_steps)) max(log10(time_steps))]); % Auto-adjust x limits

% Customize grid and axes
grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];  % Light gray grid color
ax.GridAlpha = 0.8;  % Grid transparency
ax.GridLineStyle = '-';  % Dashed grid lines
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.LineWidth = 2;
ax.FontSize = 30; % Set font size for axes

% Add legend
legend('Location', 'northwest', 'FontSize', 35, 'Interpreter', 'latex');

% Add title
title('$h = 0.02$, $T = 1$', ...
    'FontSize', 32, 'Interpreter', 'latex');

% Save figure as EPS with transparent background
exportgraphics(gcf, 'BiofilmContraction.eps', 'BackgroundColor', 'none');
