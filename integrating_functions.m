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
%% Exercise 1: Antiderivative of Polynomial
syms x C

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

%% Linearity of integration
clear; clc; close all;
syms x C

df = x^2;
% Indefinite integrals
fx1 = int(pi*df, x) + C;
fx2 = pi* int(df, x) + C;
pretty(fx1)
pretty(fx2)

% Definite integrals 
a = 0;
b = 3*pi;
fx3 = int(pi*df, x, a, b) + C;
fx4 = pi*int(df, x, a, b) + C;
pretty(fx3)
pretty(fx4)

%% Geometric Intituion
clear; clc; close all;

x = linspace(0, 5*pi/2, 500);
fx1 = x.^2;
fx2 = pi * x.^2;

figure(1); 
hold on


% Plot curves
plot(x, fx1, 'r-', 'LineWidth', 1.5);
plot(x, fx2, 'b-', 'LineWidth', 2);
legend('f(x) = x^2', 'f(x) = \pix^2', 'Location', 'northwest');
% Separate fills for each function
area(x, fx2, 'FaceColor', "b", 'FaceAlpha', 0.5, 'EdgeColor', 'none');  % Blue fill for πx2
area(x, fx1, 'FaceColor', "r", 'FaceAlpha', 0.5, 'EdgeColor', 'none');  % Pink fill for x2

% Formatting
xlim([0, 5*pi/2]);
xticks(0:pi/2:5*pi/2);
xlabel('x'); ylabel('y');

hold off


%% Summation in integration
clear; clc; close all;
syms x
fx = x^2;
gx = 10*sin(x);
hx = x^2 + 10*sin(x);
disp('Indefinite integrals')
pretty(int(fx, x) + int(gx, x))
pretty(int(hx, x))

% definite
a = 0;
b = pi;
disp('Definite integrals')
pretty(int(fx, x, a, b) + int(gx, x, a, b))
pretty(int(hx, x, a, b))

%% Visualize the summed functions
clear; clc; close all;

x = linspace(-2, 11, 500);
a = 0;
b = 10;
[~,idxa] = min(abs(x-a));
[~,idxb] = min(abs(x-b));

fx = x.^2;
gx = 10*sin(x);
hx = x.^2 + 10*sin(x);


% "trapz" calculates the numerical integral with trapezoids
int_f = trapz(x(idxa:idxb), fx(idxa:idxb));
int_g = trapz(x(idxa:idxb), gx(idxa:idxb));
int_h = trapz(x(idxa:idxb), hx(idxa:idxb));

figure(1);

subplot(1,3,1)
plot(x, fx, 'r-', 'LineWidth', 1.5);
hold on
area(x(idxa:idxb), fx(idxa:idxb), 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off
title(['\int f dx = ', num2str(int_f, '%.2f')], 'Color', 'r');

subplot(1,3,2)
plot(x, gx, 'b-', 'LineWidth', 1.5);
hold on
area(x(idxa:idxb), gx(idxa:idxb), 'FaceColor', 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
yline(0, 'b--');
hold off
title(['\int g dx = ', num2str(int_g, '%.2f')], 'Color', 'b');

subplot(1,3,3)
plot(x, hx, 'm-', 'LineWidth', 1.5);
hold on
area(x(idxa:idxb), hx(idxa:idxb), 'FaceColor', 'm', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off
title(['\int h dx = ', num2str(int_h, '%.2f')], 'Color', 'm');