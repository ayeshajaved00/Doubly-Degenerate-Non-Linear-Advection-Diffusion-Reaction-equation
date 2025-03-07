clc
clear all;
close all;
time_steps = [0.00001,0.000025,0.00005,0.0001,0.00025,0.0005,0.001, 0.0025,0.005, 0.01, 0.025,0.05,0.1,0.25,0.5,1];

%%  err <1e-10 & N=800, err = sqrt(trapz(x,(S_iter - S_intiter).^2  + dt * (dW_dx.^2)));

Arithmetic_mean_M = [6.35671e-7,1.6028e-06,3.2761e-06,7.0061e-06,2.341634723651494e-05,7.440152854770656e-05,2.769358044240997e-04,9.258538634640114e-04,0.004134199774409,0.014087518970363,0.010047043938337,0.014438305705427,0.093909631252034,0.414494867498202,0.736752599006474,0.917264123362926];

Arithmetic_mean_Adap = [6.4827e-7,1.6344e-6,3.341405544978202e-06,7.145608366864710e-06,2.387942461096346e-05,7.584393262647567e-05,2.820779261329343e-04,0.001815541550623,0.004067098617908,0.020380208101793,0.047830331323731,0.083150216070803,0.004452467430727,0.022592485712638,0.004520739124625,0.014111513037301];

Arithmetic_mean_N = [6.4831e-7,1.6348e-6,3.341406143255379e-06,7.145608897114732e-06,2.387942462735983e-05,7.584393274780392e-05,2.820779266029416e-04,0.001815541553573,0.004067098612430,0.020380208086731,0.047830331151456,0.083150215863419,0.004452466654877,0.040638514310973,0.102239771923512,0.085825498147144];

hFig = figure;
set(hFig, 'Position', [100, 100, 800, 600]); % Adjust the position and size as needed
plot(log10(time_steps), log10(Arithmetic_mean_M), 'LineStyle', '-.', 'Color', [0, 0.447, 0.741], 'LineWidth', 5, 'Marker', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.85 0.47 0.32], 'MarkerFaceColor', 'b', 'DisplayName', 'M scheme');
 hold on;
plot(log10(time_steps), log10(Arithmetic_mean_N), 'LineStyle', '-.', 'Color', [0.85, 0.325, 0.098], 'LineWidth', 5, 'Marker', '+', 'MarkerSize', 15, 'MarkerEdgeColor', [0.1 0.25 0.89], 'MarkerFaceColor', 'r', 'DisplayName', 'N scheme');
plot(log10(time_steps), log10(Arithmetic_mean_Adap), 'LineStyle', ':', 'Color', [0.466, 0.674, 0.188], 'LineWidth', 5, 'Marker', 'd', 'MarkerSize', 10, 'MarkerEdgeColor', [0.9 0.5 0.77], 'MarkerFaceColor', 'g', 'DisplayName', 'MAdap scheme');

% Set axis labels with LaTeX formatting
xlabel('$\log_{10}(\tau)$', 'FontSize', 20, 'Interpreter', 'latex');
ylabel('$\log_{10}(\alpha)$', 'FontSize', 20, 'Interpreter', 'latex');

% Add gridlines with customization
grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];  % Light gray grid color
ax.GridAlpha = 0.9;              % Grid line transparency
ax.GridLineStyle = '-';          % Solid grid lines

% Adjust axis properties
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.Box = 'on';  % Keep the box around the plot

% Set y-axis limits
ylim([-8 1]);

% Add legend
legend('Location', 'southeast', 'FontSize', 36, 'Interpreter', 'latex');

% Customize the background color
ax.Color = [1, 1, 1];

% Add title with LaTeX formatting
title('$h = 0.02$, $T = 1$', 'FontSize', 30, 'Interpreter', 'latex');

% Set font size and line width for better readability
set(gca, 'FontSize', 27);
ax.LineWidth = 1.5;

% Export figure with a transparent background
exportgraphics(gcf, 'PME_Cont1.eps', 'BackgroundColor', 'none');
