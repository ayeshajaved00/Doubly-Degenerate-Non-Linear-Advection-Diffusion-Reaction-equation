clear all
close all
 clc
 
 %% Parameter and function definitions

%van Genuchten parameter
 
global Ks
global Sp
global lambda

lambda = 0.8;

%Permeability and saturation function

Ks = @(u, lambda) (sqrt(u).*(1-(1-u.^(1/lambda)).^lambda).^2);
Bp = @(u,p, lambda) min(Ks(u,lambda)/Sp(p, lambda),1);
bp = @(u,p, lambda) min(Sp(p, lambda)/Ks(u,lambda),1);

% Defining global arrays of u-values

global du
global u_star
global u_start

du = 1e-5;
u_start=1e-5;
u_values =u_start:du:1;

% Calculate pressure for each u
p_values = 1 - ((u_values.^(-1/lambda) - 1).^(1 - lambda));

%% Find corresponding $w=\Phi(u)$-values 

w_values(1)=0;
b_values(1)=0;
B_values(1)=0;

%Integrate over the subintervals [0, p(i)]

for i = 2:length(p_values)
    % For the subsequent subintervals [p(i-1), p(i)]
    hp = (p_values(i) - p_values(i-1));
    w_values(i) = w_values(i-1) + (hp/8) *((Ks(u_values(i-1), lambda) + 3*Ks((2*u_values(i-1)+u_values(i))/3 , lambda) + 3*Ks((u_values(i-1)+2*u_values(i))/3,lambda) + Ks(u_values(i), lambda)));   
end

%Find $\Phi'(u)$ from the w-values
Phip = gradient(w_values,u_values);

%Find $u^*$ and $\Phi(u^*)$
u_star = min(find(Phip>1));
phi_ustar = w_values(u_star);
u_star = u_values(u_star);

%% calculate $b(s)$ and $B(s)$ corresponding to $s$ values
global s
global s_end
global N
%The s-input string 

s_end=2;
s = u_start:du:s_end;
N = length(s);
b_values = zeros(1,N);
B_values = zeros(1,N);

%Computing $B$-values corresponding to the s-values by using $B(s)=s$ for
%$s<u^*$ and $B(s)=\Phi(u^*)+ (s-u^*)$
for i = 1:N
    if s(i) < u_star
        B_values(i) = w_values(i);
    else
        B_values(i) = phi_ustar + s(i) - u_star;
    end
end
j = 2;
w_at_1 = w_values(end);

for i = 1:N
        while (w_values(j) < B_values(i)) && (B_values(i) < w_at_1)
            j = j + 1;
        end
        % Linearly interpolate for b(s)
        b_values(i) = (u_values(j) * (B_values(i) - w_values(j-1)) + ...
                       u_values(j-1) * (w_values(j) - B_values(i))) / ...
                      (w_values(j) - w_values(j-1));
end



%% Plotting and storing the results
%  % Create the plot for Phi(s)
%  plot(u_values, w_values,  'DisplayName', '$\Phi(s)$', 'Color', [0, 0.9, 0], 'LineWidth', 4);
%  hold on;
% % Create the plot for b(s)
%  plot(s, b_values, 'LineStyle', ':', 'Color', [0.002, 0.09876, 0.41], 'DisplayName', '$b(s)$', 'LineWidth', 4);
% 
% % % Create the plot for B(s)
%  plot(s, B_values, 'LineStyle', '-.', 'DisplayName', '$B(s)$', 'Color', [0.9, 0.09876, 0.41],'LineWidth', 4);
%  line([1 1], [w_at_1 s_end],  'Color', [0, 0.9, 0], 'LineWidth', 4);
%  xlabel('$s$', 'FontSize', 20, 'Interpreter', 'latex');
%  ylabel('Function values', 'FontSize', 20, 'Interpreter', 'latex');
% %title('$m = 6$', 'interpreter', 'latex', 'FontSize', 18);
% 
% % Display the grid
%  grid on;
% % 
% % % Customize grid appearance
% ax = gca;
% ax.GridColor = [0.7, 0.7, 0.7];  % Light gray grid color
% ax.GridAlpha = 0.5;  % Grid line transparency
% ax.GridLineStyle = '-';  % Solid grid lines
% ax.XAxisLocation = 'bottom';
% ax.YAxisLocation = 'left';
% ax.Box = 'on';  % Remove the box around the plot
% 
% % Set the background color to a light gray
% ax.Color = [1, 1, 1];
% % 
% % % Customize legend
% % % Customize legend, excluding the vertical line
%  h = legend('show', 'Location', 'Best', 'Interpreter', 'latex', 'FontSize', 26);
% % 
% % % Remove the legend entry for the vertical line
%  h.String = {'$\Phi(s)$','$b(s)$', '$B(s)$'};
% % % Increase the font size of axis labels and tick labels
% % set(gca, 'FontSize', 16);
% % 
% % % Make the axis and grid lines more prominent
%  ax.LineWidth = 1.5;
% % % Add text annotation to the plot at coordinates (x_coord, y_coord)
% % % x_coord = 1.5; % Adjust as needed
% % % y_coord = 1;  % Adjust as needed
% % % text(x_coord, y_coord, 'h', 'FontSize', 14, 'Color', 'k');
% % % Save the plot as an image with a transparent background
% %  exportgraphics(gcf, 'Richards_plot.eps', 'BackgroundColor', 'none');
% % 
