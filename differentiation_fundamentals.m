%% Limits and functions are combined into differentiation. 
% For instance, derivatives are the slopes of lines as the limit of change approaching 0.
%
% Derivative : The sequence of instantaneous slopes. (Slope of a function in each point)
% Differentiation: The process of finding the derivative (Like an algorithm/ technique to find derivative)
% 

%% Slope of a line
% f(x) = m*x + b
% m = y2 - y2 / x2 - x1 = delta y / delta x
% But we take multiple slopes as x2-x1 approaching 0 for each step of a
% segment
% 
% 1- Slope is great for a straight line
% 2- Slope for a piecewise or curved line is informative but coarse!
% 3- Compute slopes of sections (segments), each segment going smaller and
% smaller 
% 4- Derivative is the "slope series" (vector of slopes) as segment width goes to 0. 

%% Global Slopes
% Plot the (-1, 1) and (3,6) and the line in between.
% Compute the slope of that line and report the slope in the title.
clear
p1 = [-1 ,1];
p2 = [3, 6];

m = (p1(2) - p2(2)) / (p1(1) - p2(1));

% p3 = [0 , a] : intercept = a
% m = (p3(2) - p1(2)) / (p3(1) - p1(1));
% 1.2 = (a - 1) / (0 - -1);
% 1.2 = a -1
% a = 2.2 = intercept!
a = 2.2;

syms x
f = m*x + a;

figure(1);
fplot(f, [-1, 3], 'LineWidth', 2)
hold on
plot(p1(1),p1(2), 'o', 'LineWidth', 3)
plot(p2(1),p2(2), 'o', 'LineWidth', 3)
grid on
ylim([-2 7]);
xlim([-2 6]);
legend;
title("The slope of the line is :", num2str(m))
hold off
%% Local and Global Slopes
% Clean the environment
clear;
clc;
close all;

% 1. Define function and domain
N = 10;
x = linspace(-1, 5, N);
y = x.^2;

% 2. Compute the slope of each line segment
delta_y = diff(y); % Calculate the change in y (Δy)
delta_x = diff(x); % Calculate the change in x (Δx)
slopes = delta_y ./ delta_x; % Slope = Δy / Δx

% 3. Create the plots
figure(1);

% Plot the original function
subplot(2, 1, 1);
plot(x, y, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
title('Function');
xlabel('x');
ylabel('y=f(x)');
grid on;

% Plot the "slope series"
subplot(2, 1, 2);
hold on;
for i = 1:length(slopes)
    plot([x(i), x(i+1)], [slopes(i), slopes(i)], 'LineWidth', 3);
end
title('Segment slopes');
xlabel('x');
ylabel('Local slope (m)');
grid on;
ylim([-1, 9]); % Set y-axis limits to match the example
hold off;

% 4. Compare (Δy) with (Δy/Δx) in the Command Window
fprintf('--- Results for N = 5 ---\n');
disp('Δy values:');
disp(delta_y');
disp('Slope (Δy/Δx) values:');
disp(slopes');

% Compute the global slope using first and last point
delta_x_global = x(1) - x(N);
delta_y_global = y(1) - y(N);
slope_global = delta_y_global / delta_x_global;

% Compute the average of local slopes
avg_local_slopes = mean(slopes);

fprintf('Global Slope:  %s\n', num2str(slope_global));
fprintf('Average Local Slope:  %s\n', num2str(avg_local_slopes));


%% Derivative Definition 
% m = Δy / Δx 
% y = f(x)
% Δy = f(x + Δx) - f(x)
% 
% m = (f(x + Δx) - f(x)) / Δx
%   = (f(x + h) - f(x)) / h ----- We want this h (Δx) is smallest... 
% Sooo we take limit(m) h goes to 0...
% This limit is the "DERIVATIVE": = dy / dx 

% Slope when Δx is large..
% Derivative when Δx is "infinitesimal"...

