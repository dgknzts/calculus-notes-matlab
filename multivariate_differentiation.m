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

%% Gradients / Grandient Fields
% Gradient is a vector each elemant corresponds to partial derivative
% of each variable
% Nabla (reverse triangle), pronounced as "del"...
% f_x , f_y <-- partial derivatives
% Nabla = vector(f_x, f_y) <-- gradients

% + Compact representation of all partial derivatives
% + Can be used in matrix calculations (e.g., linear algeba, differential equations)
% + Used to comput directional derivatives
% + Points to the direction of the biggest increase in the function's "local"
% landscape

%% Exercise 1:
clear; clc; close all;
syms x y
fxy = x*exp(-(x^2 + y^2)); %fxy = exp(-abs(x*y));
dfdx = diff(fxy, x, 1);
dfdy = diff(fxy, y, 1);

[X, Y] = meshgrid(linspace(-3,3,100));
funcs = {fxy, dfdx, dfdy};
titles = {'f(x,y)', 'f_x', 'f_y'};

figure(1);
for i = 1:3
    subplot(1,3,i)
    Z = double(subs(funcs{i}, {x,y}, {X,Y}));
    imagesc([-3 3], [-3 3], Z)
    axis xy
    xlabel('x'); ylabel('y')
    title([titles{i} ' = ' char(funcs{i})])
end

% Gradient field!!!!!!!!!!
% Numeric grid
[X, Y] = meshgrid(linspace(-2, 2, 20));

% Convert symbolic derivatives to numeric
Zx = double(subs(dfdx, {x, y}, {X, Y}));
Zy = double(subs(dfdy, {x, y}, {X, Y}));

figure(2);
fcontour(fxy, [-2 2 -2 2], "Fill","on")
hold on
quiver(X, Y, Zx, Zy, 'black')
hold off
xlabel('x'); ylabel('y')
title(['f(x,y) = ' char(fxy)])

%% Exercise: Climb to the peak!
clear; clc; close all;
% Develop an algorithm to find the local max.
%close all; clear; clc;
syms x y
fxy = x*exp(-(x^2 + y^2));

% Numeric grid
[X, Y] = meshgrid(linspace(-3, 3, 50));
Xi = randi([1 length(X)]);
Yi = randi([1 length(Y)]);

% Num of steps
maxStep = 100;

% Convert symbolic derivatives to numeric
Zf = double(subs(fxy, {x, y}, {X, Y}));
saveXidx = zeros(length(X), 1);
saveYidx = zeros(length(Y), 1);

for i = 1:maxStep
    saveXidx(i) = Xi;
    saveYidx(i) = Yi;
    funcValue(i) = Zf(Xi, Yi);

    % Try four different directions:
    if Xi ~= length(X)
        p1 = Zf(Xi+1, Yi) - funcValue(i);
    else
        p1 = NaN;
    end

    if Xi ~= 1
        p2 = Zf(Xi-1, Yi) - funcValue(i);
    else
        p2 = NaN;
    end

    if Yi ~= length(Y)
        p3 = Zf(Xi, Yi+1) - funcValue(i);
    else
        p3 = NaN;
    end

    if Yi ~=1
        p4 = Zf(Xi, Yi-1) - funcValue(i);
    else
        p4 = NaN;
    end


    [best, idx] = max([p1, p2, p3, p4]);
    if best <= 1.0e-16
        break
    end
    
    % Update the locationnn
    if idx == 1
        Xi = Xi+1;
    elseif idx == 2
        Xi = Xi-1;
    elseif idx == 3
        Yi = Yi+1;
    elseif idx == 4
        Yi = Yi-1;
    end
    

end

figure(1);
imagesc([-3 3], [-3 3], Zf)
axis xy
hold on

colors = hot(i);
for j = 1:i
    plot(X(saveXidx(j), saveYidx(j)), Y(saveXidx(j), saveYidx(j)), 'o', ...
        'MarkerSize', 10, 'MarkerFaceColor', colors(j,:), 'MarkerEdgeColor', 'k')
end
hold off

xlabel('x'); ylabel('y')
title(sprintf('Reached peak in %d steps', i))

%% Now with using gradients
clear; clc; close all;
syms x y
fxy = x*exp(-(x^2 + y^2));
dfdx = diff(fxy, x);
dfdy = diff(fxy, y);

N = 51;
[X, Y] = meshgrid(linspace(-3, 3, N));
Xi = randi([1 N]);  % row = y direction
Yi = randi([1 N]);  % col = x direction

maxStep = 100;

Zf = double(subs(fxy, {x, y}, {X, Y}));
Zx = double(subs(dfdx, {x, y}, {X, Y}));
Zy = double(subs(dfdy, {x, y}, {X, Y}));

saveXidx = zeros(maxStep, 1);
saveYidx = zeros(maxStep, 1);

for i = 1:maxStep
    saveXidx(i) = Xi;
    saveYidx(i) = Yi;
    
    % gx controls Y index (column), gy controls X index (row)
    gx = Zx(Xi, Yi);  % ∂f/∂x → move in column direction
    gy = Zy(Xi, Yi);  % ∂f/∂y → move in row direction
    
    if abs(gx) < 1e-10 && abs(gy) < 1e-10
        break
    end
    
    % Move based on gradient
    if abs(gx) > abs(gy)
        if gx > 0, Yi = Yi + 1; else, Yi = Yi - 1; end  % x direction = column
    else
        if gy > 0, Xi = Xi + 1; else, Xi = Xi - 1; end  % y direction = row
    end
    
    % Boundary fix
    Xi = max(1, min(N, Xi));
    Yi = max(1, min(N, Yi));
    
    % Stuck check
    if i > 1 && saveXidx(i) == Xi && saveYidx(i) == Yi
        break
    end
end

% Plot
figure(1);
imagesc([-3 3], [-3 3], Zf)
axis xy
hold on

colors = parula(i);
for j = 1:i
    plot(X(saveXidx(j), saveYidx(j)), Y(saveXidx(j), saveYidx(j)), 'o', ...
        'MarkerSize', 10, 'MarkerFaceColor', colors(j,:), 'MarkerEdgeColor', 'k')
end
hold off

xlabel('x'); ylabel('y')
title(sprintf('Gradient climb: %d steps', i))
colorbar

%% m-D Gradient Descent
% Algorithm-concept-math all same
% Perform 1D gradient descent on each direction independently.
% Its convenient to store all directions in one vector instead of sep vars.
clear; clc; close all;
syms x y
fxy = 3*(1-x)^2*exp(-x^2-(y+1)^2) - 10*(x/5 - x^3 - y^5) * exp(-x^2-y^2)-11/3*exp(-(x+1)^2-y^2);

delf_x = diff(fxy, x);
delf_y = diff(fxy, y);

figure(1);
title(['f(x,y) = ' char(fxy)])
subplot(1,3,1)
fcontour(fxy, [-3 3 -3 3], 'Fill','on')
xlabel('x'); ylabel('y')
subplot(1,3,2)
fcontour(delf_x, [-3 3 -3 3], 'Fill','on')
xlabel('x'); 
subplot(1,3,3)
fcontour(delf_y, [-3 3 -3 3], 'Fill','on')
xlabel('x');

nStep = 100;
[X, Y] = meshgrid(linspace(-3, 3, nStep));
Zx = double(subs(delf_x, {x, y}, {X, Y}));
Zy = double(subs(delf_y, {x, y}, {X, Y}));

% Gradient descent
learnRate = 0.01;
nEpoch = 500;

% Random start point
px = 6*rand(1) - 3;
py = 6*rand(1) - 3;

% Store history
histX = zeros(1, nEpoch);
histY = zeros(1, nEpoch);

for i = 1:nEpoch
    histX(i) = px;
    histY(i) = py;
    
    % Get gradient at current position
    gx = double(subs(delf_x, {x, y}, {px, py}));
    gy = double(subs(delf_y, {x, y}, {px, py}));
    
    % Update position (descend = subtract gradient)
    px = px + learnRate * gx;
    py = py + learnRate * gy;
end

fprintf('Final position: (%.4f, %.4f)\n', px, py)
fprintf('Final f(x,y): %.4f\n', double(subs(fxy, {x,y}, {px,py})))

% Plot result
figure(2);
fcontour(fxy, [-3 3 -3 3], 'Fill','on')
hold on
plot(histX, histY, 'black', 'LineWidth', 3, 'MarkerSize', 5)
plot(histX(1), histY(1), 'go', 'MarkerSize', 7, 'MarkerFaceColor', 'g')  % Start
plot(histX(end), histY(end), 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r')  % End
hold off
xlabel('x'); ylabel('y')
legend('', 'Path', 'Start', 'End')
colorbar


%% Find exact solution with critical values
xCritical = solve(delf_x == 0); % WE CANT....
% It is even two variables... IRL there are thousands of variables...

%% IRL We need empirical solutions!!!!!!!!!!!!!



