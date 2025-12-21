%% Application of L'Hopital's:  Racing functions to infinitiy
% Which function grows faster to inf:
% to determine we can ratio them as
% lim x->inf f(x) / g(x) .

% Possible answers:
% inf = Numerator faster
% 0 = Denominator faster
% r = growth ratio settles to a value
% !!! Use L-Hopital's again and again to solve them!


% Note: The natural exponent grows to infinity faster than
% any power (or any non-infinitely differentiable function)


%% Second Derivative Test
% Purpose: distinguish minimum and maximum as critical points

% if ~x is a critical point of f, and
% if f'(~x) and f''(~x) exist, then:

% f''(~x) > 0 -> f(~x) is a local min
% f''(~x) < 0 -> f(~x) is a local max
% f''(~x) = 0 -> inconclusive

%% Compute and draw until second derivative
clear; clc; close all;

syms x
f = 2*pi*x + sin(pi*x) ; % we can test different functions here

secDerTest(f, x)

function secDerTest(f, x)
    df1 = diff(f, x);
    df2 = diff(f, x, 2);
    
    df1_zeros = double(solve(df1 == 0, x, 'Real', true));
    
    if isempty(df1_zeros)
        warning('No critical points found!')
        return
    end
    
    f_str = ['f = ' char(f)];
    df1_str = ['f'' = ' char(df1)];
    df2_str = ['f'''' = ' char(df2)];
    
    f_at_zeros = double(subs(f, x, df1_zeros));      % Eklendi
    df2_at_zeros = double(subs(df2, x, df1_zeros));
    
    figure(1);
    fplot(f, 'LineWidth', 2)
    hold on
    fplot(df1, 'LineWidth', 2)
    fplot(df2, 'LineWidth', 2)
    ylim([-40 40])
    grid on
    yline(0, '--k', 'LineWidth', 1)
    plot(df1_zeros, f_at_zeros, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
    plot(df1_zeros, df2_at_zeros, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', [0.5 0 0.5])
    hold off
    legend(f_str, df1_str, df2_str, 'y=0', 'f''=0', 'f''''=0')
    
    for i = 1:length(df1_zeros)
        xi = df1_zeros(i);
        df2_val = double(subs(df2, x, xi));  % double() ile çevir
        
        fprintf('Point %d: x = %.4f\n', i, xi);
        if df2_val > 0
            fprintf('  Result: Local MINIMUM\n\n');
        elseif df2_val < 0
            fprintf('  Result: Local MAXIMUM\n\n');
        else
            fprintf('  Result: INCONCLUSIVE\n\n');
        end
    end
end


%% Linear Approximations
% finding the value of the function x0 is the goal.
% find a point near the target: point a.
% compute the derivative at point a. extend that derivative to the x0.
% it is a better estimate compared to the near target a = x0...

% formula: L(x0) = f(a) + f'(a)*(x0-a)
% L is linear approximation
% intercept = f(a)
% slope = f'(a) 
% dx = x0-a (distance)

% This might fail due to the distance (dx). But it depends on the function.
% What is acceptable error is also depends on the field.
% Approximation over or under estimate? Depends on the curvature, we can do
%   second derivative test.

%% Example
clear; close all; clc;
syms x real % put real to ignore complex nums

f(x) = abs(sin(x)); %sqrt(x) + 2*x;
a = (5*pi)/8;
x0 = pi/2;

L = linApprox(f, x, a, x0);
disp(L)
analytic_x0 = double(f(x0));
disp(analytic_x0)

% % Do it for multiple a's
a = linspace(-2*pi, 2*pi, 50);
x0 = mean(a);
Ls = zeros(size(a));

for i = 1:length(a)
    Ls(1,i) = linApprox(f, x, a(i), x0);
end

errors = double(subs(f, x, x0))- Ls;

figure(1);
subplot(2,1,1)
plot(a, errors, '-o', 'LineWidth', 2);
xline(x0, '--k', 'LineWidth', 1);
subplot(2,1,2)
fplot(f, 'LineWidth', 2)
hold on
plot(a, Ls, '-o', 'LineWidth', 2)
xline(x0, '--k', 'LineWidth', 1);
hold off

function L = linApprox(f, x, a, x0)
    df = diff(f, x, 1);
    fa = double(subs(f, x, a));
    dfa = double(subs(df, x, a));
    L = dfa*(x0 - a) + fa;
end

%% Newton's Method for finding roots
% Goal: Finding the x when y = 0. 
% We start with an x0 point, and draw a tangent line and find the point in
% that line y = 0. And this will be the next point x1 we will draw a tangent.

% Bad initial guess can lead to lack of convergence (this is usuall in ML)
% - Pick a initial guess closer to f(xn) = 0
% f'(xn) = 0
% - Pick another starting point
% Unknown or difficulty derivative
% - Use the "secant method" empirical local slope
% F(x) gas no real valued roots
% - Accept this answer, or find minimum value

%% Exercise: Find the real valued roots with solve and newtons method
clear; clc; close all;
syms x
f(x) = cos(x) - x^2; % Change function forfun
actRoots = solve(f(x), x, "Real", true);
disp(double(actRoots))

% Implement Newton's Method
nIter = 10;
startPoint = pi; % If you make it negative you can find the second root..
[newtonRoot, eachIter] = newtonMethod(f, x, startPoint, nIter);
disp(newtonRoot)

figure(1);
plot(1:nIter, eachIter, '-o', 'LineWidth', 2)
yline(double(actRoots), '--r', 'LineWidth', 2);
xlabel('Iteration')
ylabel('x value')
legend('Newton iterations', 'Actual root')

figure(2);
fplot(f, 'LineWidth', 2);
ylim([-10 10])
xlim([-1 2])
yline(0, '--r', 'LineWidth', 2);
xline(double(actRoots), '--b', 'LineWidth', 2);
hold on
for p = 1:nIter
    plot(eachIter(p), double(f(eachIter(p))), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g')
end
hold off
legend('f(x)', 'y=0', 'Root', 'Iterations')

% Newton function
function [out, answers] = newtonMethod(f, x, startPoint, nIter)
    answers = zeros(1, nIter);
    for i = 1:nIter
        fx = double(subs(f, x, startPoint));
        dfx = double(subs(diff(f), x, startPoint));
        out = startPoint - fx / dfx;
        startPoint = out;
        answers(1, i) = out;
    end
end

%% Example optimizaiton problem
% One wants to fence in a rectangle of 400m^2. Shares one side with
% someone. Other person will pay 50% of the cost of that side. Fencinsing
% is 100 euro/meter. Find the lengths of the sides that minimze the fencing
% cost.

clear; clc; close all;
syms x

% a*b = 400
% 2*a + b + b/2 = cost function
% Reduce to one variable -> a = 400/b
f(x) = 800/x + x + x/2;
df = diff(f, x);
df2 = diff(f, x, 2);

min_cost = solve(df, x, "Real", true);
secondDerivativeTest = double(subs(df2, x, min_cost));

% Filter only minimums (df2 > 0)
valid_mins = min_cost(secondDerivativeTest > 0);
optimal_x = double(valid_mins);

% Calculate results
optimal_a = 400 / optimal_x;  % a = 400/b
total_cost = double(subs(f, x, optimal_x)) * 100;  % euro

% Print results
fprintf('Side b (shared): %.2f meters\n', optimal_x)
fprintf('Side a: %.2f meters\n', optimal_a)
fprintf('Area: %.2f m²\n', optimal_x * optimal_a)

% Plot
figure(1);
fplot(f, [5 150], 'LineWidth', 2)
hold on
plot(optimal_x, double(f(optimal_x)), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
hold off
xlabel('Side b (meters)')
ylabel('Cost (×100 euro)')
title('Fencing Cost Optimization')
legend('Cost function', 'Minimum')
grid on

%% Gradient Descent
% 1) Guess a solution
% 2) Compute the error (cost, minimization objective)
% 3) Learn from mistakes to improve the next guess

% Problems: cannot find the exact value, local/global minimum etc..

% Gradient descent vs f' = 0
% GD is better for high dimensional problems (e.g., millions of variables)
% Used when derivative is unknown or can only be locally estimated
% "Good enough" solutions are usually okay in ML
% Most limitations have solutions or mitigations

%% Implement GD
syms x
f(x) = 3*x^2 - 3*x + 4;
df = diff(f, x);


nIter = 100;
learningRate = .01;
startPoint = 4*rand(1);


allEpochs = zeros(1,nIter);
for i = 1:nIter
    allEpochs(1, i) = startPoint;
    direction = subs(df, x, startPoint);
    startPoint = startPoint - direction * learningRate;
end

double(startPoint)

% Plot
figure(1);
fplot(f, [-2 2], 'b', 'LineWidth', 2)
hold on
fplot(df, [-2 2], 'Color', [0.85 0.33 0.1], 'LineWidth', 2)
plot(allEpochs, double(f(allEpochs)), 'ko', 'MarkerFaceColor', 'k')
plot(allEpochs, double(subs(df, x, allEpochs)), 'ko', 'MarkerFaceColor', 'k')
hold off

xlabel('x'); 
ylabel('f(x)');
title(sprintf('Local minimum: %.4f', double(startPoint)))
legend('f(x)', 'df')