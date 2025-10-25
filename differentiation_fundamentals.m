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

%% Derivative of a line
syms x
f = (5/4)* x + 9/4;

figure(1)
figure('Color','white');
fplot(f, [-5 6]);
xlabel('x');
ylabel('f(x)');
ax = gca;
ax.XAxisLocation = 'origin'; % Center the axes at (0,0)
ax.YAxisLocation = 'origin';
ax.XTick = -5:1:5;
box off; % Remove the box outline

aa = diff(f);
fprintf('\n Derivative of f(x): %s\n', char(aa));

%% Another one
clear;
clc;
close all;

syms x
f = x^2;
bb = diff(f);
fprintf('\n Derivative of f(x): %s\n', char(bb));

% Plot the derivative
figure(1);
fplot(bb);
xlabel('x');
ylabel('f(x)');
ax = gca;
ax.XAxisLocation = 'origin'; % Center the axes at (0,0)
ax.YAxisLocation = 'origin';
ax.XTick = -5:1:5;
box off; % Remove the box outline


deriv_at_neg1 = subs(bb, x, -1);
deriv_at_0 = subs(bb, x, 0);
deriv_at_2 = subs(bb, x, 2);

% Print the results
fprintf('\n Derivative of x = -1: %s\n', char(deriv_at_neg1));
fprintf('\n Derivative of x = 0: %s\n', char(deriv_at_0));
fprintf('\n Derivative of x = 2: %s\n', char(deriv_at_2));

%% Derivative of summed terms
clear; clc; close all;

syms x a b c
f1(x) = a*x^2;
f2(x) = b*x^3;
f3(x) = c*exp(2*x );
combinedf(x) = f1(x) + f2(x) + f3(x); 

% YOU CAN ALSO SELECT COMPONENTS OF A FUNCTION WITH THIS!!!!!!
%components = children(combinedf);
%disp(components{1})
pretty(combinedf(x))

% Derivatife of combined function
j = diff(combinedf, x);
pretty(j)
% Derivative of ax^2
pretty(diff(f1(x)))
% Derivativ of b*x^3
pretty(diff(f2(x)))
% Derivative of c* exp(2*x)
pretty(diff(f3(x)))



%% Derivatives of trig functions
% sin'(x) = cos(x)
% cos'(x) = -sin(x)

clear; clc; close all;

% Create x-variable from -1.5π to 1.5π
dx = 0.1;
x = -1.5*pi:dx:1.5*pi;

% Calculate cos(x)
y = cos(x);

% Calculate the discrete derivative using diff
dy_dx = diff(y) / dx;

% Since diff reduces array size by 1, adjust x for the derivative
x_diff = x(1:end-1);

% Calculate -sin(x) for comparison
neg_sin = -sin(x);

% Create the plot
figure(1);
plot(x, y, 'LineWidth', 1, 'DisplayName', 'cos(x)');
hold on;
plot(x_diff, dy_dx, 'LineWidth', 1, 'DisplayName', 'diff(cos(x))');
plot(x, neg_sin, 'o','LineWidth', 2, 'MarkerSize', 5, 'DisplayName', '-sin(x)');
hold off;

xlabel('Angle (rad.)');
ylabel('Value');
legend('show');
grid on;
title('Exercise 1: d/dx cos(x)');

%% Do the exact same for (d/dx)* sin(x)
clear; clc; close;

dx = 0.1;
x = -2*pi:dx:2*pi;

fx = sin(x);
fx_diff = diff(fx) / dx;
x_diff = x(1:end-1);

figure(1);
plot(x, fx, 'LineWidth', 2);
grid on;
hold on;
plot(x_diff, fx_diff);
plot(x, cos(x), 'o');

