%% Limits via Zeno's Paradox
% This script demonstrates the concept of a limit by iteratively
% approaching a point from both the left and the right.

% Define the function using an anonymous
% f(x) = cos^2(x^2) + pi
f = @(x) cos(x.^2).^2 + pi;

% --- Plotting---
figure(1);
clf; % Clear the current figure
fplot(f, [-2.2, 2.2], 'LineWidth', 2);
grid on;
hold on; % Hold the plot to draw on top of it

%% Approximate the limit of the function via Zeno's paradox
% We will find the limit: lim x->1 f(x) = ?

a = 1;        % The point to approach
nRep = 10;    % Number of iterations

% Starting x-values (on either side of 'a')
x0 = [a-1, a+1]; % Starts at [0, 2]

% Pre-allocate matrices to store the results for each iteration
xAxisVals = zeros(nRep, 2);
limitVals = zeros(nRep, 2);

% Zeno's paradox algorithm
for i = 1:nRep
    % Store the current x and f(x) values
    xAxisVals(i, :) = x0;
    limitVals(i, :) = f(x0);
    
    % Update the x-values (move them halfway towards 'a')
    x0 = (x0 + a) / 2;
end

% --- Plot the approximation steps ---
% Plot the left approach points (in blue)
plot(xAxisVals(:, 1), limitVals(:, 1), 'b:o', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');

% Plot the right approach points (in green)
plot(xAxisVals(:, 2), limitVals(:, 2), 'r:o', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');

% Update the title and legend
title(['Limit Approximation: lim x->1 f(x) ≈ ' num2str(f(a))]);
legend('f(x)', 'Left Approach', 'Right Approach', 'Location', 'northwest');
hold off;

%% Another limit aproximation 
% 1. Setup
clear; clc; close all;
a = 3;      % Point to approach
nRep = 15;  % Number of steps

% 2. Define Function and Run Algorithm
f = @(x) cos(x*pi) + log(x).^2;
x_hist = zeros(nRep, 2);
y_hist = zeros(nRep, 2);
x_points = [a - 1, a + 1]; % Starting points

for i = 1:nRep
    x_hist(i, :) = x_points;
    y_hist(i, :) = f(x_points);
    x_points = (x_points + a) / 2; % Move halfway closer to 'a'
end

% 3. Visualize and Display Results
figure('Color','white');
hold on;

fplot(f, [a - 1, a + 1], 'k', 'LineWidth', 2);
plot(x_hist(:, 1), y_hist(:, 1), 's', 'MarkerFaceColor', '#0072BD', 'MarkerSize', 8); % Left
plot(x_hist(:, 2), y_hist(:, 2), 'o', 'MarkerFaceColor', '#D95319', 'MarkerSize', 8); % Right

grid on;
title(sprintf('Function value at x=%d is %.15f', a, f(a)));
legend('f(x)', 'From the left', 'From the right', 'Location', 'NorthWest');
hold off;

fprintf('Limit from left:  %.15f\n', y_hist(end, 1));
fprintf('Limit from right: %.15f\n', y_hist(end, 2));
fprintf('Actual value:     %.15f\n', f(a));

%% Find the limits with "symbolic" variable and function
syms x
f_sym = cos(x*pi) + log(x)^2;

% 2. Calculate the analytic limits
limit_left = limit(f_sym, x, 3, 'left');
limit_right = limit(f_sym, x, 3, 'right');
actual_value = subs(f_sym, x, 3); % Substitute x=3 into the function

% 3. Display the symbolic (exact) results
disp('Limit as x approaches 3 from the left:');
disp(limit_left);

disp('Limit as x approaches 3 from the right:');
disp(limit_right);

disp('Exact function value at x=3 (Symbolic):');
disp(actual_value);

% 4. Convert the symbolic result to a numerical value
numerical_value = double(actual_value);

fprintf('\nFunction value at limit (Numeric):\n%.15f\n', numerical_value);


%% Infinite Limits Example
clear; clc; close all;

%% 1. Symbolic Limit Calculation

% Define symbolic variable and function
syms x
f_sym = 1 / (x-2)^2;

% Calculate limits as x approaches 2
lim_left  = limit(f_sym, x, 2, 'left');
lim_right = limit(f_sym, x, 2, 'right');
lim_total = limit(f_sym, x, 2);

% Display results
fprintf('Limit as x->2 from the left:  %s\n', char(lim_left));
fprintf('Limit as x->2 from the right: %s\n', char(lim_right));
fprintf('Limit as x->2 (two-sided):    %s\n', char(lim_total));


%% 2. Plotting the Function

figure('Color', 'white');
hold on;

% Plot the function itself
fplot(f_sym, [0, 4], 'LineWidth', 2, 'Color', '#0072BD');

% Draw the vertical asymptote at x=2
xline(2, '--', 'Color', [0.5 0.5 0.5]);

% Set plot limits and labels
ylim([0 100]);
xlim([0 4]);
grid on;
box on;
title('$f(x) = \frac{1}{(x-2)^2}$', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('x');
ylabel('f(x)');

hold off;

%% Undefined Limits Example
clear; clc; close all;

%% 1. Symbolic Limit Calculation

% Define the symbolic variable and function
syms x
f_sym = abs(x - 2) / (x - 2);

% Calculate the one-sided limits
lim_left  = limit(f_sym, x, 2, 'left');
lim_right = limit(f_sym, x, 2, 'right');
lim_total = limit(f_sym, x, 2); % Will be NaN (undefined)

% Display the results
fprintf('Limit as x->2 from the left:  %s\n', char(lim_left));
fprintf('Limit as x->2 from the right: %s\n', char(lim_right));
fprintf('Limit as x->2 (two-sided):    %s\n', char(lim_total));


%% 2. Plotting the Function

figure('Color','white');
% fplot is excellent at finding and showing discontinuities
fplot(f_sym, [-2, 4], 'LineWidth', 2);

% --- Customize the plot to look like the example ---
ax = gca;
ax.XAxisLocation = 'origin'; % Center the axes at (0,0)
ax.YAxisLocation = 'origin';
box off; % Remove the box outline
ylim([-1.5, 1.5]);
xlabel('x');
ylabel('f(x)');
title('$f(x) = \frac{|x-2|}{x-2}$', 'Interpreter', 'latex', 'FontSize', 16);

%% Properties of limits
% Factorable: lim(cf) = c*lim(f)
% Additive: lim(f + g) = lim(f) + lim(g)
% Multiplicative: lim(f*g) = lim(f) * lim(g)
% Divisive: lim(f / g) = lim(f) / lim(g)

%% Exercise: Limits are factorable
clear; clc; close all;
% Define the function
syms x
f_x = x^3/3 + 100*sqrt(abs(x));

% Show that as x -> 5, lim(cf) = c(lim(f)) 
% Define random constant c
c = rand(1);

lim_left_1  = limit(c*f_x, x, 5, 'left');
lim_right_1 = limit(c*f_x, x, 5, 'right');
lim_total_1 = limit(c*f_x, x, 5); % Will be NaN (undefined)

lim_left_2  = c*limit(f_x, x, 5, 'left');
lim_right_2 = c*limit(f_x, x, 5, 'right');
lim_total_2 = c*limit(f_x, x, 5); % Will be NaN (undefined)

disp(double(lim_left_1) == double(lim_left_2));
disp(double(lim_right_1) == double(lim_right_2));
disp(double(lim_total_1) == double(lim_total_2));

%% Exercise: Limits are additive
clear; clc; close all;
syms x

fx = log(x) + x^2;
gx = exp(-x) + x^3;


lim_total_1 = limit(fx, x, pi) + limit(gx, x, pi);
lim_total_2 = limit(fx+gx, x, pi);

disp(double(lim_total_1));
disp(double(lim_total_2));

%% Exercise: Limits are multiplicative
lim_total_1 = limit(fx, x, pi) * limit(gx, x, pi);
lim_total_2 = limit(fx*gx, x, pi);
disp(double(lim_total_1));
disp(double(lim_total_2));

%% More:
lim_total_1 = limit(fx, x, pi) * limit(fx, x, pi) * limit(fx, x, pi);
lim_total_2 = limit(fx^3, x, pi);
disp(double(lim_total_1));
disp(double(lim_total_2));

%% Exercise: Limits are divisive
clear; clc; close all;
syms x

fx = log(x) + x^2;
gx = exp(-x) + x^3;
hx = x^3 + x^2 + x;

%% as x -> pi:
lim_1 = limit(fx, x, pi) / limit(gx, x, pi);
lim_2 = limit(fx/gx, x, pi);
disp(double(lim_1));
disp(double(lim_2));

%% as x -> 0:
lim_3 = limit(gx, x, 0) / limit(hx, x, 0);
lim_4 = limit(gx/hx, x, 0);
disp(double(lim_3));
disp(double(lim_4));
%% The problem here: These properties are valid only when the limit exist!! 
% There is no limit for two sided limit in the lim4....
lim_3 = limit(gx, x, 0) / limit(hx, x, 0, "right");
lim_4 = limit(gx/hx, x, 0, "right");
disp(double(lim_3));
disp(double(lim_4));

lim_3 = limit(gx, x, 0) / limit(hx, x, 0, "left");
lim_4 = limit(gx/hx, x, 0, "left");
disp(double(lim_3));
disp(double(lim_4));
% Plot
figure('Color','white');
% fplot is excellent at finding and showing discontinuities
fplot(gx/hx, [-5, 5], 'LineWidth', 2);

% --- Customize the plot to look like the example ---
ax = gca;
ax.XAxisLocation = 'origin'; % Center the axes at (0,0)
ax.YAxisLocation = 'origin';
box off; % Remove the box outline
ylim([-1.5, 1.5]);
xlabel('x');
ylabel('f(x)');

%% Defining Discontinuities with limits.
% Jump Discontinuity: One sided limits exist but differ
% Infinite Discontinuity: One or two-sided limits are Inf.
% Removable Discontinuity: Limit is not same with the func value

%% Jump Discontinuity
% Estimate the limits for the piecewise function
clear; clc; close all;

% 1. Define the three pieces of the function IN MATLAB FUNCTION!!
function y = piecewise_f(x)
%   The function is defined as:
%     - y = sin(x*pi)      if x < 0
%     - y = 1.5            if x = 0
%     - y = -(x - 2)^2     if x > 0

% Define a small tolerance to handle floating-point comparisons
tol = 1e-6; % A small number for the "bullseye" around zero

% Initialize the output vector
y = zeros(size(x));

% Create logical indices based on the tolerance
idx_neg = (x < -tol);
idx_pos = (x > tol);
idx_zero = (abs(x) < tol); % Check if x is "close enough" to zero

% Apply the function's rules to the corresponding elements
y(idx_neg) = sin(x(idx_neg) * pi);
y(idx_zero) = 1.5;
y(idx_pos) = -(x(idx_pos) - 2).^2;

end

% --- Example Script Using piecewise_f ---
% 1. Define the domain
x_values = linspace(-1, 2, 1000);

% 2. Call your new function to get the y-values
y_values = piecewise_f(x_values);

% 3. Create the plot
figure('Color', 'white');
plot(x_values, y_values, 'o', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', 'none');

grid on;
xlabel('x');
ylabel('y=f(x)');
title('Plot Generated from piecewise_f.m');


% --- Limit Estimation Algorithm ---
nRep = 10; % Number of estimation steps

% For the left limit (approaching 0 from negative side)
x_left = -0.1; % Starting point
x_hist_left = zeros(nRep, 1);
y_hist_left = zeros(nRep, 1);
for i = 1:nRep
    x_hist_left(i) = x_left;
    % The function automatically uses the correct rule because x_left < 0
    y_hist_left(i) = piecewise_f(x_left);
    x_left = x_left / 2; % Move halfway to 0
end

% For the right limit (approaching 0 from positive side)
x_right = 0.1; % Starting point
x_hist_right = zeros(nRep, 1);
y_hist_right = zeros(nRep, 1);
for i = 1:nRep
    x_hist_right(i) = x_right;
    % The function automatically uses the correct rule because x_right > 0
    y_hist_right(i) = piecewise_f(x_right);
    x_right = x_right / 2; % Move halfway to 0
end

% --- Display Text Results ---
fprintf('Limit approaches %.15f from the left.\n', y_hist_left(end));
fprintf('Limit approaches %.15f from the right.\n', y_hist_right(end));
fprintf('Function value at x=0: %.1f\n\n', piecewise_f(0));

% --- Create the Plot ---
figure('Color','white');
hold on;
plot(x_hist_left, y_hist_left, 's', 'MarkerFaceColor', '#0072BD', 'MarkerSize', 8);
plot(x_hist_right, y_hist_right, 'o', 'MarkerFaceColor', '#D95319', 'MarkerSize', 8);
% Get the actual value directly from the function
plot(0, piecewise_f(0), 'rx', 'MarkerSize', 12, 'LineWidth', 2);

% Customize plot
grid on;
title(sprintf('Function value at x=0 is %.1f.', piecewise_f(0)));
xlabel('x');
ylabel('y');
legend('From the left', 'From the right', 'y=f(0)', 'Location', 'NorthEast');
hold off;

%% Jump Discontinuity
%% Symbolic analysis of a jump discontinuity
clear; clc; close all;

% 1. Define the symbolic variable
syms x

% 2. Create the symbolic piecewise function
f_sym = piecewise( ...
    x < 0, sin(x*pi), ...
    x == 0, 1.5, ...
    x > 0, -(x-2)^2);

% 3. Calculate the limits
lim_left  = limit(f_sym, x, 0, 'left');
lim_right = limit(f_sym, x, 0, 'right');
lim_total = limit(f_sym, x, 0); % Will be NaN (undefined)

% 4. Find the actual function value by substitution
val_at_0 = subs(f_sym, x, 0);

% 5. Display all results
fprintf('Limit as x approaches 0 from the left:  %s\n', char(lim_left));
fprintf('Limit as x approaches 0 from the right: %s\n', char(lim_right));
fprintf('Two-sided limit as x approaches 0:    %s\n', char(lim_total));
fprintf('Function value at limit:              %s\n', char(val_at_0));

%% Symbolic analysis of an infinite discontinuity
clear; clc; close all;

% 1. Define the symbolic variable and function
syms x
f_sym = 3 / (1 - x^2);

% 2. Calculate limits as x approaches -1
lim_left  = limit(f_sym, x, -1, 'left');
lim_right = limit(f_sym, x, -1, 'right');
lim_total = limit(f_sym, x, -1); % Will be NaN

% 3. Display the results
fprintf('Limit as x approaches -1 from the left:  %s\n', char(lim_left));
fprintf('Limit as x approaches -1 from the right: %s\n', char(lim_right));
fprintf('Two-sided limit as x approaches -1:    %s\n', char(lim_total));

% 4. Plot the function and its asymptotes
figure('Color', 'white');
hold on;

fplot(f_sym, [-2.5, 2.5], 'LineWidth', 2);
xline(-1, '--k'); % Asymptote at x=-1
xline(1, '--k');  % Asymptote at x=1

% Customize the plot
ylim([-10, 10]);
grid on;
box on;
xlabel('x');
ylabel('f(x)');
title('$f(x) = \frac{3}{1-x^2}$', 'Interpreter', 'latex', 'FontSize', 16);

hold off;

%% Symbolic analysis of an oscillating discontinuity
clear; clc; close all;

% 1. Define the symbolic variable and function
syms x
f_sym = sin(1/x);

% 2. Calculate limits as x approaches 0
% The limit does not exist because of the infinite oscillation.
% MATLAB will return NaN (Not a Number).
lim_left  = limit(f_sym, x, 0, 'left');
lim_right = limit(f_sym, x, 0, 'right');
lim_total = limit(f_sym, x, 0);

% 3. Display the results
fprintf('Limit as x approaches 0 from the left:  %s\n', char(lim_left));
fprintf('Limit as x approaches 0 from the right: %s\n', char(lim_right));
fprintf('Two-sided limit as x approaches 0:    %s\n', char(lim_total));

% 4. Plot the function to visualize the oscillation
figure('Color', 'white');
fplot(f_sym, [-0.1, 0.1], 'LineWidth', 1.5);

% Customize the plot
ylim([-1.5, 1.5]);
grid on;
box on;
xlabel('x');
ylabel('f(x)');
title('$f(x) = \sin(1/x)$', 'Interpreter', 'latex', 'FontSize', 16);