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