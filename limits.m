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