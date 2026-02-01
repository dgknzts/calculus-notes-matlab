%% Area between two curves
clear; clc; close all;
syms x
fx = log(x);
dfx = int(fx, x);

xx = 0.02:0.02:5;
fx_l = double(subs(fx, x, xx));
dfx_l = double(subs(dfx, x, xx));
hx = fx_l - dfx_l;

figure(1);

% Find intersection points
[pks, locs] = findpeaks(-abs(hx));

% Extract region
x_region = xx(locs(1):locs(2));
y1 = fx_l(locs(1):locs(2));
y2 = dfx_l(locs(1):locs(2));
patch([x_region x_region(end:-1:1)], [y1 y2(end:-1:1)], 'yellow');
hold on;

plot(xx, fx_l, 'LineWidth', 2);
plot(xx, dfx_l, 'LineWidth', 2);
plot(xx(locs), fx_l(locs), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

grid on;
hold off;

% Area
area = int(fx - dfx, x, xx(locs(1)), xx(locs(2)));
area = double(area);

%% Parametric Curves
clear; clc; close all;

syms t

% Define the shared part
shared = exp(cos(t)) - sin(t/12)^5 - 2*cos(4*t);

% Define x(t) and y(t)
xt = shared * sin(t);
yt = shared * cos(t);

% Compute derivatives
dxt = diff(xt, t);
dyt = diff(yt, t);

% Convert to numeric (lambdify)
t_vals = linspace(0, 2*pi, 1000);
x_vals = double(subs(xt, t, t_vals));
y_vals = double(subs(yt, t, t_vals));
dx_vals = double(subs(dxt, t, t_vals));
dy_vals = double(subs(dyt, t, t_vals));


figure(1);
subplot(2,2,1);
plot(x_vals, y_vals, 'k', 'LineWidth', 1.5);
xlabel('x'); ylabel('y');
title('Parametric functions');
grid on;
axis equal;

subplot(2,2,2);
plot(t_vals, x_vals, 'LineWidth', 1.5);
hold on;
plot(t_vals, y_vals, 'LineWidth', 1.5);
xlabel('t'); ylabel('x or y');
title('Functions by t');
legend('x(t)', 'y(t)');
grid on;

subplot(2,2,3);
plot(dx_vals, dy_vals, 'k', 'LineWidth', 1.5);
xlabel('dx'); ylabel('dy');
title('Their derivatives');
grid on;
axis equal;

subplot(2,2,4);
plot(t_vals, dx_vals, 'LineWidth', 1.5);
hold on;
plot(t_vals, dy_vals, 'LineWidth', 1.5);
xlabel('t'); ylabel('dx or dy');
title('Derivatives by t');
legend("x'(t)", "y'(t)");
grid on;

%% Exercise: Length by function sampling
clear; clc; close all;

syms t

% Define parametric functions
shared = exp(cos(t)) - sin(t/12)^5 - 2*cos(4*t);
xt = shared * sin(t);
yt = shared * cos(t);

% Derivatives
dxt = diff(xt, t);
dyt = diff(yt, t);

% Convert to function handles
xt_func = matlabFunction(xt, 'Vars', t);
yt_func = matlabFunction(yt, 'Vars', t);
dxt_func = matlabFunction(dxt, 'Vars', t);
dyt_func = matlabFunction(dyt, 'Vars', t);
integrand_func = @(t) sqrt(dxt_func(t).^2 + dyt_func(t).^2);

% True length using integral (spi.quad)
length_quad = integral(integrand_func, 0, 2*pi);

% Loop over different N values
N_vals = 5:15:1000;
lengths_numpy = zeros(size(N_vals));
lengths_simpson = zeros(size(N_vals));

for i = 1:length(N_vals)
    N = N_vals(i);
    t_vals = linspace(0, 2*pi, N);
    x_vals = xt_func(t_vals);
    y_vals = yt_func(t_vals);

    % Method 1: Numpy (sum of segment lengths)
    dx = diff(x_vals);
    dy = diff(y_vals);
    lengths_numpy(i) = sum(sqrt(dx.^2 + dy.^2));

    % Method 2: Simpson's rule
    f_vals = integrand_func(t_vals);
    n = N - 1;
    if mod(n, 2) ~= 0       % Simpson needs even intervals
        f_vals = f_vals(1:end-1);
        n = n - 1;
    end
    h = 2*pi / n;
    lengths_simpson(i) = h/3 * (f_vals(1) + f_vals(end) + ...
        4*sum(f_vals(2:2:end-1)) + 2*sum(f_vals(3:2:end-2)));
end

% Plotting
figure(1);

% Left: Example curves for N=10 and N=1000
subplot(1,2,1);
t_small = linspace(0, 2*pi, 10);
t_large = linspace(0, 2*pi, 1000);
plot(xt_func(t_small), yt_func(t_small), 'o-', 'LineWidth', 1.5, 'MarkerSize', 4);
hold on;
plot(xt_func(t_large), yt_func(t_large), 'LineWidth', 1.5);
xlabel('x'); ylabel('y');
title('Example curves for two Ns');
legend('N = 10', 'N = 1000');
axis equal; grid on;

% Right: Curve length convergence by N
subplot(1,2,2);
plot(N_vals, lengths_numpy, 'b-s', 'MarkerSize', 3, 'LineWidth', 1.5);
hold on;
plot(N_vals, lengths_simpson, 'x-', 'Color', [1 0.5 0], 'MarkerSize', 3, 'LineWidth', 1.5);
yline(length_quad, 'k--', 'LineWidth', 1.5);
xlabel('N'); ylabel('Curve length');
title('Curve length by N');
legend('Numpy', "Simpson's", 'spi.quad');
grid on;
ylim([35, 36])

%% The function and its solid
clear; clc; close all;

syms x
fx = log(x);
gx = int(fx, x);

disp('f(x) = '); disp(fx)
disp('g(x) = '); disp(gx)

%% Find intersection points
fx_l = matlabFunction(fx);
gx_l = matlabFunction(gx);

xx = linspace(0.01, 4.5, 523);
fun_diff = abs(fx_l(xx) - gx_l(xx));
[~, peekz] = findpeaks(-fun_diff);
intersection_points = xx(peekz);

% Draw
figure(1);
plot(xx, fx_l(xx), 'LineWidth', 2); hold on;
plot(xx, gx_l(xx), 'LineWidth', 2);
x4area = linspace(intersection_points(1), intersection_points(2), 79);
fill([x4area fliplr(x4area)], [fx_l(x4area) fliplr(gx_l(x4area))], 'k', 'FaceAlpha', 0.5);
plot(intersection_points, fx_l(intersection_points), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'w');
xlabel('x'); ylabel('y');
legend('f(x) = log(x)', 'g(x) = xlog(x) - x');
grid on;

%% Shift functions up
min2add = round(min(gx_l(xx)));
gx = int(fx, x) - min2add;
fx = log(x) - min2add;

fx_l = matlabFunction(fx);
gx_l = matlabFunction(gx);

disp('Shifted f(x) = '); disp(fx)
disp('Shifted g(x) = '); disp(gx)

%% Draw shifted
figure(2);
plot(xx, fx_l(xx), 'LineWidth', 2); hold on;
plot(xx, gx_l(xx), 'LineWidth', 2);
fill([x4area fliplr(x4area)], [fx_l(x4area) fliplr(gx_l(x4area))], 'k', 'FaceAlpha', 0.5);
plot(intersection_points, fx_l(intersection_points), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'w');
xlabel('x'); ylabel('y');
ylim([-0.1 3]);
grid on;

%% 3D Solid of Revolution
a = intersection_points(1);
b = intersection_points(2);
xx = linspace(a, b, 200);

[X, Theta] = meshgrid(xx, linspace(0, 2*pi, 100));
Y_f = fx_l(X) .* cos(Theta);
Z_f = fx_l(X) .* sin(Theta);
Y_g = gx_l(X) .* cos(Theta);
Z_g = gx_l(X) .* sin(Theta);

figure(3);
subplot(1,2,1);
plot(xx, fx_l(xx), 'LineWidth', 2); hold on;
plot(xx, gx_l(xx), '--', 'LineWidth', 2);
fill([xx fliplr(xx)], [fx_l(xx) fliplr(gx_l(xx))], 'm', 'FaceAlpha', 0.2);
xlabel('x'); ylabel('y = f(x) or g(x)');
title('2D area to revolve');
legend('f(x)', 'g(x)', 'Area to revolve');

subplot(1,2,2);
surf(X, Y_f, Z_f, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
surf(X, Y_g, Z_g, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Solid in 3D');
axis equal;

%% Exercise 2: Approximate volume and surface area
dfx_l = matlabFunction(diff(fx, x));
dgx_l = matlabFunction(diff(gx, x));

dx = xx(2) - xx(1);

% Volume via Simpson's rule
volume = pi * simpson(fx_l(xx).^2 - gx_l(xx).^2, dx);

% Surface areas
surfacearea_f = 2*pi * simpson(fx_l(xx) .* sqrt(1 + dfx_l(xx).^2), dx);
surfacearea_g = 2*pi * simpson(gx_l(xx) .* sqrt(1 + dgx_l(xx).^2), dx);

fprintf('Volume: %.2f (au)^3\n\n', volume);
fprintf('Inner surface area: %.2f (au)^2\n', surfacearea_g);
fprintf('Outer surface area: %.2f (au)^2\n', surfacearea_f);
fprintf('Total surface area: %.2f (au)^2\n', surfacearea_f + surfacearea_g);

%% Simpson's rule function (like scipy.integrate.simpson)
function I = simpson(y, dx)
    n = length(y) - 1;
    if mod(n, 2) ~= 0
        y = y(1:end-1);
        n = n - 1;
    end
    I = dx/3 * (y(1) + y(end) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2)));
end