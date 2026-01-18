%% Improper Integrals
clear; clc; close all;
syms x

% Define functions and limits
f1 = 1/x^2;      a1 = 1;    b1 = inf;
f2 = pi/x^3;     a2 = 1;    b2 = inf;
f3 = x^4;        a3 = -inf; b3 = -1;
f4 = 1/sqrt(x);  a4 = 1;    b4 = inf;

% Calculate integrals
I1 = int(f1, x, a1, b1);
I2 = int(f2, x, a2, b2);
I3 = int(f3, x, a3, b3);
I4 = int(f4, x, a4, b4);

%% Exercise 2
clear; clc; close all;
syms x

% Define range of p values
p_values = -2:0.2:2;

% Create figure with two subplots
figure(1)
% Left plot: Indefinite integral functions
subplot(1,2,1)
x_plot = linspace(1, 10, 200);
hold on
for p = p_values
    % Compute indefinite integral (antiderivative)
    f = x^-p;
    F = int(f, x);
    
    % Convert to numeric function and plot
    F_numeric = matlabFunction(F);
    y_plot = F_numeric(x_plot);
    
    % Color based on convergence (p < -1 converges)
    if p < -1
        plot(x_plot, y_plot, 'g', 'LineWidth', 2)
    else
        plot(x_plot, y_plot, 'r', 'LineWidth', 2)
    end
end

xlabel('x')
ylabel('y = F(x)')
title('Indefinite integral functions')
ylim([-15 30])
grid on
hold off

% Right plot: Convergence indicator
subplot(1,2,2)
hold on
for i = 1:length(p_values)
    p = p_values(i);
    % Check convergence using definite integral
    f = x^-p;
    I = int(f, x, 1, inf);
    
    % Plot convergence status
    if isfinite(double(I))
        plot(p, 1, 'gs', 'MarkerSize', 7, 'MarkerFaceColor', 'g')  % Converges
    else
        plot(p, 2, 'rs', 'MarkerSize', 7, 'MarkerFaceColor', 'r')  % Diverges
    end
end

xlabel('p in x^(-p)')
title('Integral convergence')
yticks([1 2])
yticklabels({'Converges', 'Diverges'})
ylim([0.5 2.5])
xlim([min(p_values)-0.5, max(p_values)+0.5])
grid on
hold off