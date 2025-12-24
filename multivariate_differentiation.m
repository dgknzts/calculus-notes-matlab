%% Partial Derivatives

% What? 
% Multivariable functions change along each variable. Partial derivatives
% reveal those changes

% Why?
% We need to decompose function changes to understand the functions'
% bheaviors.

% How?
% Differentiate one variable at a time, treating the other variables like
% constants. 

% The other variable exist in the partial deriative of other variable is
% says that the change in one depends on the other variable... 
% Effect on changing x depends on the y.

% Formal definition:
% fx = delf / delx = lim h ->0 f(x + h, y) - f(x, y) / h
% h is only added to x. y is act like a constant.

%% Exercise 1: Differentiate a multivariable function
syms x y
fxy = x^2 / pi^3 + sin(y);

dfdx = diff(fxy, x);
pretty(dfdx);

dfdy = diff(fxy, y);
pretty(dfdy);

% Visualize % fmesh, fsurf, fcontour all cool!!!!!!!
figure(1);
subplot(3,1,1)
[X, Y] = meshgrid(linspace(0,6,100));
Z = double(subs(fxy, {x,y}, {X,Y})); % Convert to matrix with a domain
imagesc([0 6], [0 6], Z)
colorbar
colormap('jet')
axis xy
xlabel('x'); ylabel('y'); zlabel('f(x,y)')
title('f(x,y) = x^2/\pi^3 + sin(y)')
subplot(3,1,2)
[X, Y] = meshgrid(linspace(0,6,100));
Z = double(subs(dfdx, {x,y}, {X,Y})); % Convert to matrix with a domain
imagesc([0 6], [0 6], Z)
colorbar
colormap('jet')
axis xy
xlabel('x'); ylabel('y'); zlabel('f(x,y)')
title('fx = 2x/pi^3')
subplot(3,1,3)
[X, Y] = meshgrid(linspace(0,6,100));
Z = double(subs(dfdy, {x,y}, {X,Y})); % Convert to matrix with a domain
imagesc([0 6], [0 6], Z)
colorbar
colormap('jet')
axis xy
xlabel('x'); ylabel('y'); zlabel('f(x,y)')
title('fx = cos(y)')

%% Exercise
clear; clc; close all;

syms x y 
fxy = x^2*sin(y) + x^2 + sin(y);

dfdx = diff(fxy, x);
pretty(dfdx);

dfdy = diff(fxy, y);
pretty(dfdy);

figure(1);
subplot(3,1,1)
[X, Y] = meshgrid(linspace(0,6,100));
Z = double(subs(fxy, {x,y}, {X,Y})); % Convert to matrix with a domain
imagesc([0 6], [0 6], Z)
colorbar
colormap('jet')
axis xy
xlabel('x'); ylabel('y'); zlabel('f(x,y)')
title('f(x,y) = x^2*sin(y) + x^2 + sin(y)')
subplot(3,1,2)
[X, Y] = meshgrid(linspace(0,6,100));
Z = double(subs(dfdx, {x,y}, {X,Y})); % Convert to matrix with a domain
imagesc([0 6], [0 6], Z)
colorbar
colormap('jet')
axis xy
xlabel('x'); ylabel('y'); zlabel('f(x,y)')
title('f_x = 2*x + 2*x*sin(y)')
subplot(3,1,3)
[X, Y] = meshgrid(linspace(0,6,100));
Z = double(subs(dfdy, {x,y}, {X,Y})); % Convert to matrix with a domain
imagesc([0 6], [0 6], Z)
colorbar
colormap('jet')
axis xy
xlabel('x'); ylabel('y'); zlabel('f(x,y)')
title('f_y = cos(y)*x^2 + cos(y)')

%% Higher-order partial derivatives
clear; clc; close all;

% Notatin: del* (delf / delx) / delx --> Second order
% Or: del^2f / del*x^2 
% Or: (f_x)_x
% Or : f_xx % And if the second der is y -> f_xy

% Is f_xy = f_yx .
% If first derivatives are continous at all x,y points
% it is symmetric. (Clairaut's Theorem)
syms x y
fxy = cos(2*x*y^2) + x^2 - y*exp(x);

df_x = diff(fxy, x);
disp(df_x)

df_y = diff(fxy, y);
disp(df_y)

df_xx = diff(fxy, x, 2);
disp(df_xx)

df_yy = diff(fxy, y, 2);
disp(df_yy)

df_xy = diff(df_x, y, 1);
disp(df_xy)

df_yx = diff(df_y, x, 1);
disp(df_yx)

%% Exercise: Random Function Generator
clear; clc; close all;
syms x y

% Random coefficient (-5 to 5, excluding 0)
coefPool = [-5:-1, 1:5];
randCoef = @() coefPool(randi(length(coefPool)));

% Random power (-5 to 5)
randPow = @() randi([-5, 5]);

% Transcendental functions
transFuncs = {@cos, @sin, @log, @exp};
randTrans = @() transFuncs{randi(4)};

% Term 1: Polynomial in x * y
c1 = randCoef();
p1 = randPow();
term1 = c1 * x^p1 * y;

% Term 2: Polynomial in y * x
c2 = randCoef();
p2 = randPow();
term2 = c2 * y^p2 * x;

% Term 3: Transcendental with xy
c3 = randCoef();
transFunc = randTrans();
term3 = c3 * transFunc(x*y);

% Final function
fxy = term1 + term2 + term3;

% Display
fprintf('\nf(x,y) = %s\n', char(fxy))

% Plot
figure;
fsurf(fxy, [-2 2 -2 2])
xlabel('x'); ylabel('y'); zlabel('f(x,y)')
title(['f(x,y) = ' char(fxy)])
colorbar

disp('Second order partial derivatives:')
disp(['df_xx = ' char(diff(fxy, x, 2))])
disp(['df_yy = ' char(diff(fxy, y, 2))])