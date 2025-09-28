%% MATLAB Functions Basics - Comprehensive Mathematical Function Analysis
% Converted from functions_basics.mlx Live Script
% This script demonstrates various mathematical functions used in calculus
% with professional visualization techniques and educational examples

clear; clc; close all;

%% Section 1: Basic Polynomial Function
% Mathematical Function: y = x² + 3x³ - x⁴
% Demonstrates element-wise operations and basic plotting

x = -2:0.01:2;  % Domain with high resolution (0.01 step)
y = x.^2 + 3*x.^3 - x.^4;  % Polynomial function with element-wise operations

figure(1);
plot(x, y, 'LineWidth', 2);
title('Basic Polynomial Function: $y = x^2 + 3x^3 - x^4$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 12);
legend('$y = x^2 + 3x^3 - x^4$', 'Interpreter', 'latex', 'FontSize', 12);
grid on;

%% Section 2: Coordinate System Visualization
% Mathematical Function: f(β) = -β⁴ + 3β³ + β² (same function, different notation)
% Demonstrates custom axis positioning with origin at center

beta = -2:0.01:2;
f_beta = -beta.^4 + 3*beta.^3 + beta.^2;

figure(2);
plot(beta, f_beta, 'LineWidth', 2, 'Color', 'red');
title('Mathematical Coordinate System: $f(\beta) = -\beta^4 + 3\beta^3 + \beta^2$', ...
       'Interpreter', 'latex', 'FontSize', 14);
xlabel('$\beta$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$f(\beta)$', 'Interpreter', 'latex', 'FontSize', 12);

% Center the coordinate system at origin
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
grid on;

%% Section 3: Random Polynomial Generation
% Mathematical Function: f(x) = a₀ + a₁x + a₂x² + a₃x³
% Shows dynamic coefficient generation and formatted titles

x = -4:0.1:4;  % Extended domain for third-order polynomial

% Generate random coefficients
a0 = rand(1) * 10;
a1 = rand(1) * 10;
a2 = rand(1) * 10;
a3 = rand(1) * 10;

% Calculate polynomial with random coefficients
f_random = a0 + a1*x + a2*x.^2 + a3*x.^3;

figure(3);
plot(x, f_random, 'LineWidth', 2, 'Color', 'green');
title(sprintf('Random Third-Order Polynomial: f(x) = %.2f + %.2fx + %.2fx² + %.2fx³', ...
              a0, a1, a2, a3), 'FontSize', 12);
xlabel('x', 'FontSize', 12);
ylabel('f(x)', 'FontSize', 12);
grid on;

%% Section 4: Rational Functions
% Mathematical Function: f(x) = (x² - 2x)/(x² - 4)
% Demonstrates asymptotic behavior and high-resolution plotting

x = linspace(-5, 5, 1000);  % High resolution for smooth asymptotes

% Rational function with vertical asymptotes at x = ±2
numerator = x.^2 - 2*x;
denominator = x.^2 - 4;
f_rational = numerator ./ denominator;

figure(4);
plot(x, f_rational, 'LineWidth', 2, 'Color', 'blue');
title('Rational Function: $f(x) = \frac{x^2 - 2x}{x^2 - 4}$', ...
       'Interpreter', 'latex', 'FontSize', 14);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 12);

% Center coordinate system and set limits
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
ylim([-10, 10]);
grid on;

%% Section 5: Maclaurin Series Approximation
% Mathematical Concept: Maclaurin series for sin(x) = Σ((-1)^(n+1) * x^(2n-1)) / (2n-1)!
% Demonstrates series convergence and function approximation

x = -2*pi:0.1:2*pi;
y_exact = sin(x);  % Exact sine function

% Calculate 10-term Maclaurin series approximation
y_approx = zeros(size(x));
for n = 1:10
    term = ((-1)^(n+1) .* x.^(2*n-1)) / factorial(2*n-1);
    y_approx = y_approx + term;
end

figure(5);
plot(x, y_exact, 'LineWidth', 2, 'Color', 'red');
hold on;
plot(x, y_approx, '--', 'LineWidth', 2, 'Color', 'blue');
hold off;

title('Maclaurin Series Approximation of sin(x)', 'FontSize', 14);
xlabel('x', 'FontSize', 12);
ylabel('y', 'FontSize', 12);
legend('Exact: sin(x)', '10-term Maclaurin Series', 'Location', 'best');
grid on;

%% Section 6: Exponential Function Analysis
% Purpose: Estimates Euler's number (e) and analyzes convergence
% Mathematical Formula: lim(n→∞) (1 + 1/n)^n = e
% Demonstrates convergence analysis and error calculation

n_values = 1:100;  % Number of iterations
y_convergence = zeros(size(n_values));

% Calculate convergence to e
for i = 1:length(n_values)
    n = n_values(i);
    y_convergence(i) = (1 + 1/n)^n;
end

% Error analysis
e_exact = exp(1);
error = e_exact - y_convergence;

figure(6);
subplot(2,1,1);
plot(n_values, y_convergence, 'LineWidth', 2, 'Color', 'magenta');
hold on;
yline(e_exact, '--', 'LineWidth', 2, 'Color', 'red');
hold off;
title('Convergence to Euler''s Number: $(1 + \frac{1}{n})^n \to e$', ...
       'Interpreter', 'latex', 'FontSize', 14);
xlabel('n', 'FontSize', 12);
ylabel('$(1 + \frac{1}{n})^n$', 'Interpreter', 'latex', 'FontSize', 12);
legend('$(1 + \frac{1}{n})^n$', '$e = 2.718...$', 'Interpreter', 'latex');
grid on;

subplot(2,1,2);
plot(n_values, error, 'LineWidth', 2, 'Color', 'black');
title('Error: $e - (1 + \frac{1}{n})^n$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('n', 'FontSize', 12);
ylabel('Error', 'FontSize', 12);
grid on;

%% Section 7: Exponential Function Comparison
% Purpose: Compares different exponential function variations
% Mathematical Functions: e^(x²), e^((-x)²), e^(-x²), (e^(-x))², e^x
% Demonstrates matrix-based function storage and comprehensive comparison

x = -1:0.01:1;
nFunc = 5;  % Number of exponential functions
fVector = zeros(nFunc, length(x));

% Define different exponential functions
fVector(1, :) = exp(x.^2);        % Gaussian growth
fVector(2, :) = exp((-x).^2);     % Same as above ((-x)² = x²)
fVector(3, :) = exp(-x.^2);       % Gaussian bell curve
fVector(4, :) = (exp(-x)).^2;     % Exponential decay squared
fVector(5, :) = exp(x);           % Standard exponential

figure(7);
plot(x, fVector, 'LineWidth', 2);
title('Exponential Function Variations', 'FontSize', 14);
xlabel('x', 'FontSize', 12);
ylabel('f(x)', 'FontSize', 12);
legend({'$e^{x^2}$ (Gaussian Growth)', '$e^{(-x)^2}$ (Same as above)', ...
        '$e^{-x^2}$ (Gaussian Bell)', '$(e^{-x})^2$ (Decay Squared)', ...
        '$e^x$ (Standard Exponential)'}, ...
       'Interpreter', 'latex', 'Location', 'best');
grid on;

%% Section 8: Logarithmic vs Exponential
% Purpose: Beautiful comparison of log and exp functions
% Mathematical Functions: log(x), exp(x), log(exp(x)), exp(log(x))
% Demonstrates identity functions and careful domain handling

x_pos = 0.1:0.01:3;  % Positive domain for log(x)
x_all = -1:0.01:3;   % Extended domain for exp(x)

figure(8);
set(gcf, 'Position', [100, 100, 1000, 700]);  % Custom figure size

% Plot log(x) - only for positive x
plot(x_pos, log(x_pos), 'LineWidth', 3, 'Color', [0.8, 0.2, 0.2]);
hold on;

% Plot exp(x) - for all x
plot(x_all, exp(x_all), 'LineWidth', 3, 'Color', [0.2, 0.2, 0.8]);

% Plot identity functions
plot(x_pos, log(exp(x_pos)), '--', 'LineWidth', 2, 'Color', [0.2, 0.8, 0.2]);
plot(x_pos, exp(log(x_pos)), ':', 'LineWidth', 2, 'Color', [0.8, 0.2, 0.8]);

% Plot y = x reference line
plot(x_all, x_all, 'k--', 'LineWidth', 1);

hold off;

title('Logarithmic vs Exponential Functions', 'FontSize', 16);
xlabel('x', 'FontSize', 14);
ylabel('y', 'FontSize', 14);
legend({'$\log(x)$', '$e^x$', '$\log(e^x) = x$', '$e^{\log(x)} = x$', '$y = x$'}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northwest');
grid on;
xlim([-0.5, 3]);
ylim([-0.5, 3]);

%% Section 9: Trigonometric Functions - Basic
% Purpose: Explores basic trigonometric variations
% Mathematical Functions: sin(x), sin²(x), sin(x²)
% Demonstrates matrix-based function storage for trigonometric functions

x = -pi:0.1:2*pi;
nFunc = 3;
fVector = zeros(nFunc, length(x));

% Define trigonometric function variations
fVector(1, :) = sin(x);           % Basic sine
fVector(2, :) = sin(x).^2;        % Sine squared
fVector(3, :) = sin(x.^2);        % Sine of x squared

figure(9);
plot(x, fVector, 'LineWidth', 2);
title('Trigonometric Function Variations', 'FontSize', 14);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 12);
legend({'$\sin(x)$', '$\sin^2(x)$', '$\sin(x^2)$'}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
grid on;

%% Section 10: Advanced Trigonometric Functions
% Purpose: Demonstrates composite trigonometric functions
% Mathematical Functions: sin(cos(x)), cos(sin(x)), cos(x)
% Shows function composition and complex mathematical expressions

x = -2*pi:0.1:2*pi;
nFunc = 3;
fVector = zeros(nFunc, length(x));

% Define composite trigonometric functions
fVector(1, :) = sin(cos(x));      % Sine of cosine
fVector(2, :) = cos(sin(x));      % Cosine of sine
fVector(3, :) = cos(x);           % Basic cosine for reference

figure(10);
plot(x, fVector, 'LineWidth', 2);
title('Composite Trigonometric Functions', 'FontSize', 14);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 12);
legend({'$\sin(\cos(x))$', '$\cos(\sin(x))$', '$\cos(x)$'}, ...
       'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
grid on;

%% Section 10: Tangent Function
% Basic tangent function plotting with different resolutions

% Plot 1: Different resolutions
figure(10);
x1 = 0:0.01:2*pi;        % Coarse resolution
x2 = 0:0.05:2*pi;       % Fine resolution
y1 = tan(x1);
y2 = tan(x2);


plot(x1, y1, 'r-', 'LineWidth', 3, 'DisplayName', 'Coarse (0.05)');
hold on;
plot(x2, y2, 'b-', 'LineWidth', 2, 'DisplayName', 'Fine (0.01)');
hold off;

title('Tangent Function - Different Resolutions');
xlabel('x');
ylabel('tan(x)');
legend('Location', 'northeast');
ylim([-10, 10]);
grid on;

% Plot 2: Zoomed view
figure(11);
x_zoom = 1.55:0.001:1.60;
y_zoom = tan(x_zoom);

plot(x_zoom, y_zoom, 'LineWidth', 2, 'Color', 'red');
title('Zoomed View: tan(x) near π/2');
xlabel('x');
ylabel('tan(x)');
grid on;

%% Section 11: Piecewise Function
% Simple piecewise function plotting and evaluation

% Define domain
x = -2:0.01:5;

% Calculate piecewise function values
y = zeros(size(x));
for i = 1:length(x)
    if x(i) < 0
        y(i) = 0;
    elseif x(i) < 3
        y(i) = -2 * x(i);
    else
        y(i) = 0.1 * x(i)^3;
    end
end

% Plot
figure(12);
plot(x, y, 'LineWidth', 2);
title('Piecewise Function');
xlabel('x');
ylabel('f(x)');
grid on;

% Find specific value like your example
xloc = 0.5;
% Find the index where x is closest to xloc
xidx = find(abs(x - xloc) == min(abs(x - xloc)));
result = y(xidx);

% Display result
fprintf('xloc = %.1f\n', xloc);
fprintf('f(%.1f) = %.1f\n', xloc, result);

% Mark the point on the plot
hold on;
plot(xloc, result, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
hold off;

%% Section 12: Discontinuous Functions
% A function f(x) is continuous at x=a if the following three conditions are met;
% otherwise, the function is discontinuous at point a.

% 1. f(a) is defined.
% 2. The limit of f(x) as x approaches a exists.
% 3. The limit of f(x) as x approaches a equals f(a).


%% Jump Disconinuouty
% Define the domain D.
x = -1:0.01:2;
y = NaN(size(x));

% Loop through each x value to calculate the corresponding y.
for i = 1:length(x)
    if x(i) < 0
        y(i) = sin(x(i)*pi);
    elseif x(i) == 0
        y(i) = 1.5;
    else % x(i) > 0
        y(i) = -(x(i)-2)^2;
    end
end

% Plot the function.
figure;
plot(x(x<0), y(x<0), 'LineWidth', 2);
hold on
plot(x(x==0), y(x==0), 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');
plot(x(x>0), y(x>0), 'LineWidth', 2);

% To add the title, labels, and the single point at (0, 1.5).
title('A function with a jump discontinuity');
xlabel('x');
ylabel('y=f(x)');
hold off; % Resets the plot behavior.

%% Removable Discontinuity
% Domain
x = -1:0.1:2;
y = NaN(size(x));

% Loop through each x value to calculate the corresponding y.
for i = 1:length(x)
    if x(i) == 0
        y(i) = pi;
    else % x(i) > 0
        y(i) = sin(x(i)*pi) + x(i)^2;
    end
end

% Plot the function.
figure;
plot(x, y, 'o');


%% Infinite Discontinuity
% f(x) = 3 / (1-x^2)
% D : -2 <= x <= 2

x = -2:0.1:2;

y = 3 ./ (1-x.^2);

plot(x, y, 'LineWidth', 2)
xlim([-2,2])
ylim([-14, 15.5])


%% Solving with Symbolic
% 1. Define the symbolic variable and the entire function
syms x;
f(x) = 3 / (1 - x^2);

% 2. Automatically extract the denominator from the function
[~, denominator] = numden(f); % The ~ ignores the numerator part

% 3. Solve for where the denominator equals zero
singularities = solve(denominator == 0, x);
disp(singularities);


%% Oscillating Discontinuity
% 1. Create a numerical vector for x
x_vals = linspace(-1, 2, 200);

% 2. Calculate y
y_vals = sin(0.1 ./ (1 - x_vals));

% 3. Manually insert a break at the singularity to prevent a connecting line
y_vals(abs(x_vals - 1) < 0.01) = NaN; 

% 4. Use the standard plot function
plot(x_vals, y_vals, 'LineWidth', 1);
grid on;
xlabel('x');
ylabel('y');


%% Composite Functions
x = linspace(-5,5,100);

fx = 2.*x.^2 - 4;
gx = 7.*abs(x) + 3;

fgx = 2.*gx.^2 - 4;
gfx = 7.*abs(fx) + 3;

figure;
plot(x, fx, "LineWidth", 2);
hold on;
plot(x, gx, "LineWidth", 2);
plot(x, fgx, "LineWidth", 2);
plot(x, gfx, "LineWidth", 2);
ylim([-10,50]);
xlim([-3 3])
legend('$f(x) = 2x^2 - 4$', ...
       '$g(x) = 7|x| + 3$', ...
       '$f(g(x))$', ...
       '$g(f(x))$', ...
       'Interpreter', 'latex', ...
       'Location', 'best');
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
grid on;
hold off;

%% More about composite functions in """symbolic""" package. Better.
% 1. Define 'x' as a symbolic variable
clear;
clc;
syms x

% 2. Define f(x) and g(x) as symbolic functions
f(x) = log(2*x);
g(x) = exp(x) / 2;

% 3. Compose the functions symbolically. The syntax is very intuitive.
fgx(x) = f(g(x)); % This computes f(g(x))
gfx(x) = g(f(x)); % This computes g(f(x))

% Display the resulting symbolic expressions in the command window
fprintf('f(g(x)) = \n');
disp(fgx);

fprintf('g(f(x)) = \n');
disp(gfx);

% 4. Plot all four functions using a single fplot command
figure;
fplot([f(x), g(x), fgx(x)], [-1, 5], 'LineWidth', 2);
hold on;
fplot(gfx(x), [-1, 5], 'o', 'LineWidth', 2);
ylim([-4, 10]);
grid on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
legend('f(x)', 'g(x)', 'f(g(x))', 'g(f(x))', 'Location', 'northwest');
hold off;

%% Inverse Functions
syms x

% 2. Define f(x) and g(x) as symbolic functions
f(x) = sin(x);
g(x) = log(x);
h(x) = 2*x^2 + 5;

% 3. Compose the functions symbolically. The syntax is very intuitive.
fghx(x) = f(g(h(x))); % This computes f(g(x))

% Display the resulting symbolic expressions in the command window
fprintf('f(g(h(x))) = \n');
disp(fghx);

% 4. Plot all four functions using a single fplot command
figure;
fplot(fghx(x), [-100, 100], 'LineWidth', 2);
hold on;
ylim([-1.2, 1.2]);
grid on;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
box off;
title(sprintf('f(g(h(x))) = %s', string(fghx(x))));
hold off;

%% Computing inverse with symbolic
syms x

% 2. Define f(x) and g(x) as symbolic functions
f(x) = 2*x + 3;
f_inv(x) = finverse(f);
disp(f_inv);

f_invf(x) = f_inv(f(x));
ff_inv(x) = f(f_inv(x));

g(x) = 2*x + sin(x);
g_inv(x) = finverse(g);
disp(g_inv)

% Demonstrate that: f_inv(f(4)) = f(f_inv(4)) = 4
disp(f_inv(f(4)) == f(f_inv(4)))

figure;
fplot(f_invf(x), [-1, 10], 'LineWidth', 2);
hold on;
fplot(ff_inv(x), [-1, 10], 'o', 'LineWidth', 2 );
grid on;
xline(4, '--r', 'LineWidth', 1.5, 'Label', 'x = 4');
yline(4, '--r', 'LineWidth', 1.5, 'Label', 'y = 4');
xlabel('x');
ylabel('y');
title('Verifying f_{inv}(f(x)) = f(f_{inv}(x)) = x');
legend('f_{inv}(f(x))', 'f(f_{inv}(x))', 'Intersection at (4,4)', 'Location', 'best');
axis equal; % Ensure x and y axes have the same scale
hold off;
