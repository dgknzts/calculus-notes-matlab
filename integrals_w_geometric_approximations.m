%% Riemann Sum
clear; clc; close all;
a = -0.5;
b = pi;
x = linspace(a, b, 15);
dx = x(2) - x(1);  % Fix 1: doğru dx hesaplama
fx = cos(x) + exp(x) / 10;

figure(1);
hold on

% Draw Riemann bars (right rule)
for i = 2:length(x) 
    rectangle('Position', [x(i)-dx, 0, dx, fx(i)], 'FaceColor', 'b', 'EdgeColor', 'k')
end

syms xa
fxa = cos(xa) + exp(xa) / 10;
fplot(fxa, [-1 3.5], 'k', 'LineWidth', 2)
xlim([-1, 3.5]);
xline(a, '--g', 'LineWidth', 1.5)
xline(b, '--r', 'LineWidth', 1.5)
yline(0, '--', 'Color', [0.5 0.5 0.5])
hold off

xlabel('x')
ylabel('y = f(x)')
title(sprintf('Riemann Sum, Δx = %.2f', dx))

% Different areas
area_left = sum(fx(1:end-1)) * dx;
fprintf('Area Left ≈ %.4f\n', area_left)

area_right = sum(fx(2:end)) * dx;
fprintf('Area Right ≈ %.4f\n', area_right)

x_mid = (x(1:end-1) + x(2:end)) / 2;
fx_mid = cos(x_mid) + exp(x_mid) / 10;
area_mid = sum(fx_mid) * dx;
fprintf('Area Mid ≈ %.4f\n', area_mid)

area_anal = double(int(fxa, xa, a, b));
fprintf('Area Anal ≈ %.4f\n', area_anal)


%% Show converge with shrinking deltaX
clear; clc; close all;

a = -0.5;
b = pi;
f = @(x) cos(x) - exp(x)/10;

nBins = 5:200;
area_L = zeros(size(nBins));
area_R = zeros(size(nBins));
area_M = zeros(size(nBins));

for i = 1:length(nBins)
    x = linspace(a, b, nBins(i));
    dx = x(2) - x(1);
    
    area_L(i) = sum(f(x(1:end-1))) * dx;
    area_R(i) = sum(f(x(2:end))) * dx;
    area_M(i) = sum(f((x(1:end-1)+x(2:end))/2)) * dx;
end

plot(nBins, area_L, nBins, area_R, nBins, area_M, '--k', 'LineWidth', 2)
yline(2.7367, '--k', 'LineWidth', 2)
legend('Left', 'Right', 'Midpoint', 'Analytic')
xlabel('Number of bins')
ylabel('Area (a.u.)')
