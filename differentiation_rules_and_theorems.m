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

%% 