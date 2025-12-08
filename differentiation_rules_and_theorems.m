%% Linearity of differentiation
% Why it is linear? Proof:
% k(x) = a*f(x) + b*g(x) 

% k' = lim h -> 0 ( k(x+h) - k(x) ) / h
%    = lim h -> 0 ( (a*f(x+h) + b*g(x+h)) - a*f(x) - b*g(x) ) / h
%    = lim h -> 0 ( a*(f(x+h) - f(x)) + b*(g(x+h) - g(x)) ) / h
% k' = ( a* f' ) + ( b* g')

%% THEOREM:  Differentiability implies continuity 
% If derivative exist at point x0,
% then the function is continuous at x0.

% We start with definition of continuity: lim x->a f(x) = f(a)

% lim x->a (f(x) - f(a)) = 0  % we can do it because there is no x in f(a)
% its a constant.

% Then we multiply and divide by x-a..
% lim x->a (f(x) - f(a)) / (x - a) * lim x->a (x - a) = 0
% First part is derivative of f and second part is x = a since x -> a and 0
% f'(x) * 0 = 0 : CORRECT !

% OPPOSITE WAY IS NOT CORRECT! IF FUNCTION IS CONTINOUS IT DOES NOT MEAN IT
% HAS TO HAVE DERIVATIVE like f(x) = |x|.

%% PRODUCT RULE
% ( f(x)*g(x) )' = f'(x)*g(x) + f(x)*g'(x) 


%% CHAIN RULE
% f(g(x))' = f'(g(x))*g'(x)

% y = g(x)
% z  = f(x)
% dz/dx = dz/dy * dy/dx

%% QUOTIENT RULE
% (f/g)' = (f'g - fg') / g^2
% low d-high minus high d-low, square the bottom and away we go!

%% Confirming the product rule
syms x 
f(x) = x^2;
g(x) = cos(x);

fg = diff(f(x) * g(x));

df = diff(f(x));
dg = diff(g(x));
fg_manual = df*g(x) + f(x)*dg ;

figure(1);
fplot(fg, 'LineWidth', 2, 'Marker', 'o');
hold on
fplot(fg_manual, 'LineWidth', 2)
hold off

pretty(fg)
pretty(fg_manual)

%% Confirm the quatient rule
clear;
clc;

syms x
f(x) = x^2;
g(x) = cos(x);

fg = diff(f / g);

fg_manual = (diff(f)*g - f*diff(g)) / g^2;

pretty(fg)
pretty(fg_manual)

%% Chain rule with three functions
clear;
clc;

syms x
f(x) = 3*x^2;
g(x) = x^3 + log(x);
h(x) = cos(x);

fgh = f(g(h(x)));

fgh_diff = diff(fgh);


figure(2);
fplot(fgh, 'LineWidth', 2);
hold on
fplot(fgh_diff, 'LineWidth', 2);
xlim([-20 20]);  % Set x-axis limits
ylim([-100 100]);
hold off

%% Implicit Differentitation
% Implicit functions are the functions we cannot define without the
% circular inference of the function
% y = 2*x --> this is normal explicit.
% BUT:
% x*y = 17*x^cos(y) or sqrt(x*y) = log(sin(x*y)) + 4y^2.
% We cannnot solve for y, make y (function of x) alone on one side. 

% In implicit differentiation, we take y as y'. Then solve for y'. 
clear; clc;
syms x y

% Implicit equation xy = 1
eq = x*y - 1;

% Implicit differentiation: dy/dx = -df/dx / df/dy
y_prime = -diff(eq, x) / diff(eq, y);

% Solve for explicit form
y_explicit = solve(eq == 0, y);

% Substitute to get y' without y
y_prime_explicit = simplify(subs(y_prime, y, y_explicit));


% Plot
figure;
subplot(1,2,1)
fimplicit(eq, [-5 5 -5 5], 'b', 'LineWidth', 2)
title('xy = 1')
axis equal; grid on

subplot(1,2,2)
fplot(y_prime_explicit, [-5 5], 'r', 'LineWidth', 2)
ylim([-10 2])
title('y'' = -1/x^2')
grid on

%% Another implicit example:
clear; clc;
syms x y
eq = exp(x^2 + y^2) - x - y;

y_prime = -diff(eq, x) / diff(eq, y);
y_prime = simplify(y_prime);
pretty(y_prime)


% Plot the implicit function
figure(1);
fimplicit(eq, [-2 2 -2 2], 'b', 'LineWidth', 2)
title('e^{x^2+y^2} = x + y')
grid on; axis equal

y_explicit = solve(eq == 0, y); %Cannot solve..

%% Differentiate c^x and x^x ...
clear; clc;
syms x 
f = 2^x;
z = exp(x);
c = x^x;
f_d = diff(f);
c_d = simplify(diff(c));

pretty(f_d)
pretty(c_d)
figure(1);
fplot(f_d)
hold on
fplot(f);
fplot(z);
fplot(c);
fplot(c_d)
xlim([-1, 3])
ylim([0, 8])
title('f(x) = 2^x')
grid on; axis equal
legend;
hold off

%% HIGHER ORDER DERIVATIVES
% Second derivative: derivative of derivative
% Nth order derivative: taking the derivative N times.

% second deriv: d^2y / dx^2 in leibniz notation. or (d^2 / (dx)^2) * y
% or f'' in lagrange. for third or more: f^(3): parantheses is important.

% Tip: polynomials tend to simplify in higher order derivative.
% Tip2: transcendental functions tend to get more complicated with more diff

% Some applications: 
% curve-sketching (inflection points)
% classifying critical points
% physics (position, velocity, acceleration, jerk..)
% engineering...

%% Exercise for derivatives of derivatives
clear; clc;

syms x
f = 0.1*x^4 + exp(-x^2) + cos(x);

df = diff(f);
ddf = diff(df);
dddf = diff(ddf); % or diff(f, 3), second arg is order.

figure(1);
subplot(2,2,1)
fplot(f, "LineWidth", 2);
title("f")
ylim([-20 20])
subplot(2,2,2)
fplot(df, "LineWidth", 2);
title("f'")
ylim([-20 20])
subplot(2,2,3)
fplot(ddf, "LineWidth", 2);
title("f''")
ylim([-20 20])
subplot(2,2,4)
fplot(dddf, "LineWidth", 2);
title("f^{(3)}")
ylim([-20 20])

%% L'Hopital's Rule
% If the limit of a ratio of two functions is indeterminate, try the limit
% of the ratio of the derivatives

% If you cannot solve a limit problem --> it ended up as 0 / 0 or inf / inf
% You can take derivative of both upper and lower parts of the function
% until it is solveable...

% Assumptions before apply:
% 1) Is expressed as a ratio of functions
% 2) lim(f/g) is indeterminate
% 3) f' and g' exist
% 4) g' is not equal to 0

%% Rolle's Theorem
% if f(a) = f(c)
% then this exist: a < b < c s.t. f'(b) = 0

% f(x) must be continuous in [a, c]
% f(x) must be differentiable in (a, c)
% f(a) = f(c)

%% Mean Value Theorem
% f(b) = (f(c) - f(a)) / (c - a)    ,    a < b < c
% or =
% yc - ya / xc - xa

% f(x) must be continuous in [a, c]
% f(x) must be differentiable in (a, c)

% WHY THE NAME IS MEAN VALUE THEOREM (MVT)?
% EACH OF THIS ANSWERS SAYING THE SAME THING!:
% -- The slope over an interval equals the average of the instanteous slopes
% over that interval. 
% -- The global slope equals the average of the local slopes
% -- The slope of a secant line is the average of the slopes of all tangent
% lines. 

%% MVT Algorithm
clear; clc;
syms x

f = 2*x^2 - 3*x + 1;
a = -1;
c = 2;

[b, slope] = solveMVT(f, x, a, c);
disp(b)
disp(slope)


% Visualize the above MVT
figure(1);
fplot(f, [-2 3], 'b', 'LineWidth', 2)
hold on

% Slope line
secant = slope*(x - a) + subs(f, x, a);
tangent = slope*(x - b) + subs(f, x, b);
fplot(secant, 'g--', 'LineWidth', 2)
fplot(tangent, 'r--', 'LineWidth', 2)
plot(b, double(subs(f, x, b)), 'bo', 'MarkerSize', 10, 'LineWidth', 2)
plot(a, double(subs(f, x, a)), 'ro', 'MarkerSize', 10, 'LineWidth', 2)
plot(c, double(subs(f, x, c)), 'ro', 'MarkerSize', 10, 'LineWidth', 2)
ylim([-1, 10]);
hold off
legend('f(x) = x³', 'secant', 'tangent')
grid on


% Function solveMVT:
function [b, slope] = solveMVT(f, var, a, c)
    % Derivative of f
    df = diff(f, var);
    
    % MVT: f'(b) = (f(c) - f(a)) / (c - a)
    % Solve for b
    slope = (subs(f, var, c) - subs(f, var, a)) / (c - a);
    
    % Solve
    b_all = solve(df == slope, var);
    
    % Convert to numeric
    b_all = double(b_all);
    
    % Filter: only keep b where a < b < c
    b = b_all(b_all > a & b_all < c);
end

