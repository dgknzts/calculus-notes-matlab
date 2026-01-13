%% Initial Value Problem
clear; clc; close all; 
% Given f'(x) = ... , f(0) = .... ; find f(x).
syms x

% dfx = 2*x + 3;
% cue = 2; % When x is 1.
% fknown = 1; % initial value x

% Other example
dfx = cos(x) * x^(5/2);
cue = pi; % When x is 1.
fknown = 1; % initial value x

fx_noC = int(dfx);

f1_noC = subs(fx_noC, x, fknown);
C = cue - f1_noC;
fx = fx_noC + C;


figure(1);
fplot(fx, [1 12], 'g', 'LineWidth', 3)
hold on
offsets = linspace(-2, 2, 10);
for i = 1:length(offsets)
    fplot(fx + offsets(i), 'r', 'LineWidth', 0.5) 
end 
hold off

%% INFINITE PRACTICE PROBLEMS
clear; clc; close all;

syms x C

%% Exercise 1: Antiderivative of Polynomial

% Create random polynomial
nTerms = randi([2, 4]);
coeffs = randi([-10, 10], 1, nTerms);

% Build polynomial
expr1 = 0;
for i = 1:nTerms
    expr1 = expr1 + coeffs(i)*x^(i-1);
end

% Display problem
fprintf('Find the antiderivative:\n');
pretty(expr1);

% Solution
ans1 = int(expr1, x) + C;
fprintf('\nAnswer:\n');
%pretty(ans1);

%% Exercise 2: Definite Integral

% Create random polynomial
nTerms = randi([2, 4]);
coeffs = randi([-10, 10], 1, nTerms);

% Build polynomial
expr2 = 0;
for i = 1:nTerms
    expr2 = expr2 + coeffs(i)*x^(i-1);
end

% Create bounds
bounds = sort(randi([-3, 3], 1, 2));
a = bounds(1);
b = bounds(2);
if a == b
    b = b + 1;
end

% Display problem
fprintf('\nCompute the definite integral from %d to %d:\n', a, b);
pretty(expr2);

% Solution
ans2 = int(expr2, x, a, b);
%fprintf('\nAnswer: %s\n', char(ans2));

%% Exercise 3: Initial Value Problem

% Create random differential equation f'(x)
df = randi([-5, 5])*x + randi([-5, 5]);

% Create initial conditions
x0 = randi([-5, 5]);
fx0 = randi([-5, 5]);

% Display problem
fprintf('\nGiven f''(x) = ');
disp(df);
fprintf('and f(%d) = %d\n', x0, fx0);
fprintf('Find f(x)\n');

% Solution
fx = int(df, x) + C;
constant = solve(subs(fx, x, x0) - fx0, C);
final_answer = subs(fx, C, constant);

fprintf('\nAnswer:\n');
%pretty(final_answer);