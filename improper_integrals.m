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

%% Improper Trigonometric Integrals
clear; clc; close all;

f = @(x) cos(x.^2);  % Function handle
a = pi/4; 
blist = a:0.1:4*pi;  % Start from a

% Calculate areas using numerical integration
net_areas = zeros(1, length(blist));
tot_areas = zeros(1, length(blist));

for j = 1:length(blist)
    b = blist(j);
    net_areas(j) = integral(f, a, b);         
    tot_areas(j) = integral(@(x) abs(f(x)), a, b);  
end

% Plotting
figure(1)
subplot(2,1,1)
fplot(f, [a, blist(end)], 'LineWidth', 2)
xlabel('x')
ylabel('y = cosx^2')
grid on

subplot(2,1,2)
hold on
plot(blist, tot_areas, 'r', 'LineWidth', 2)
plot(blist, net_areas, 'b', 'LineWidth', 2)
legend('Total Area |f(x)|', 'Net Area f(x)', 'Location', 'best')
xlabel('Upper Boundary (b)')
ylabel('Area')
grid on
hold off

%% Exercise
clear; clc; close all;
syms x

% Define integrals
I1 = int(abs(cos(x)), x, 0, inf) ;
I2 = int(cos(x), x, 0, inf);
I3 = int(cos(x^2), x, 0, inf);

% Display results
fprintf('int |cos(x)| dx from 0 to ∞ = %s\n', char(I1))
fprintf('int cos(x) dx from 0 to ∞ = %s\n', char(I2))
fprintf('int cos(x²) dx from 0 to ∞ = %s\n', char(I3))

%% Empirical convergence towards infinity
clear; clc; close all;

% Parameters
p_values = linspace(-3, -1.5, 20);  % Range of p values
zeta_max = 1e5;                      % Upper bound limit

figure(1)
subplot(1,2,1)
hold on

x_plot = linspace(1, 10, 200);

for i = 1:length(p_values)
    p = p_values(i);
    y = x_plot.^p;
    
    % Color gradient from blue to red
    color = [i/length(p_values), 0, 1 - i/length(p_values)];
    plot(x_plot, y, 'Color', color, 'LineWidth', 1.5)
end

xlabel('x')
ylabel('y = f(x)')
title('Function curves')
ylim([0 1])
grid on
hold off


subplot(1,2,2)
hold on
% Range of upper bounds in log scale
zeta_values = logspace(0, 5, 100); 

for i = 1:length(p_values)
    p = p_values(i);
    
    % Compute definite integral from 1 to zeta
    % Integral of x^p = x^(p+1)/(p+1)
    areas = zeros(1, length(zeta_values));
    
    for j = 1:length(zeta_values)
        zeta = zeta_values(j);
        % Analytical solution: [x^(p+1)/(p+1)] from 1 to zeta
        areas(j) = (zeta^(p+1) - 1) / (p+1);
    end
    
    % Color gradient
    color = [i/length(p_values), 0, 1 - i/length(p_values)];
    semilogx(zeta_values, areas, 'Color', color, 'LineWidth', 1.5)
end
set(gca, 'XScale', 'log')
xlabel('Upper bound')
ylabel('Definite integral')
title('Area under f(x) = x^p')
grid on
hold off
