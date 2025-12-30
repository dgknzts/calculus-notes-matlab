%% Integration as "reverse differentitation"
% Aprox formula:
% F(x_k) = i:k sum f(x_i)*deltaX
% Integration is sum while differentation is substract

%% Exercise 1: Compute and plot integral with functions
dx = 0.01;
x = -1:dx:4;  % Domain
fx = x.^3  + 4;

[df, idf] = findDerandInt(x, dx, fx, 1);


function [outDer, outInt] = findDerandInt(x, dx, fx, n)       
    outDer = diff(fx) / dx;
    outIntTemp = cumsum(outDer) * dx;
    % Normalize integral
    [~, zeroIdx] = min(abs(outIntTemp)); 
    outIntTemp = outIntTemp - outIntTemp(zeroIdx);
    outInt = outIntTemp + fx(zeroIdx);


    figure(n);
    subplot(1,3,1)
    plot(x(1:end-1), fx(1:end-1), 'o')
    xlim([x(1), x(end)])
    subplot(1,3,2)
    plot(x(1:end-1), outDer, 'o')
    xlim([x(1), x(end)])
    subplot(1,3,3)
    plot(x(1:end-1), outInt, 'o')
    hold on
    plot(x(1:end-1), fx(1:end-1), 'red', 'LineWidth', 2)
    legend('Aprox','Orgi Func');
    hold off
    xlim([x(1), x(end)])
end

%% Integral is area
% if f(x) = c, and a and b C f(x)
% int a->b f(x)dx = Area = c* (b-a)
% But if integral is indefinite (dont have a and b)
% integral is not a one area answer it is a function for
% each value x. 
% Exercise: Geometric Approximation
clear; clc; close all;
f = @(x) x.^2 - 0.5;
a = -0.5;
b = 1;

figure(1);
% Δx = 0.2
subplot(1,2,1)
dx = 0.1;
x = a:dx:b-dx;
area = sum(f(x) * dx);

fplot(f, [a b], 'k', 'LineWidth', 2)
hold on
for i = 1:length(x)
    rectangle('Position', [x(i), min(0,f(x(i))), dx, abs(f(x(i)))], ...
        'FaceColor', [0.8 0.7 0.9], 'EdgeColor', 'k')
    plot(x(i), f(x(i)), 's', 'MarkerFaceColor', [0.5 0.3 0.7], 'MarkerSize', 8)
end
fplot(f, [a b], 'k', 'LineWidth', 2)
hold off
title(sprintf('Area = %.3f when Δx=%.2f', area, dx))
xlabel('x'); ylabel('y = x^2 - 0.5')

% Δx = 0.05
subplot(1,2,2)
dx = 0.01;
x = a:dx:b-dx;
area = sum(f(x) * dx);


hold on
for i = 1:length(x)
    rectangle('Position', [x(i), min(0,f(x(i))), dx, abs(f(x(i)))], ...
        'FaceColor', [0.8 0.7 0.9], 'EdgeColor', 'k')
    plot(x(i), f(x(i)), 's', 'MarkerFaceColor', [0.5 0.3 0.7], 'MarkerSize', 5)
end
fplot(f, [a b], 'k', 'LineWidth', 2)
hold off
title(sprintf('Area = %.3f when Δx=%.2f', area, dx))
xlabel('x'); ylabel('y = x^2 - 0.5')

%% Exercise 2: A range of Δx's
clear; clc; close all;

f = @(x) x.^2 - 0.5;
a = -1;
b = 1;

% 20 logarithmic steps from 0.5 to 0.001
dxVals = logspace(log10(0.5), log10(0.001), 20);
areas = zeros(size(dxVals));

for i = 1:length(dxVals)
    dx = dxVals(i);
    x = a:dx:b-dx;
    areas(i) = sum(f(x) * dx);
end

% Plot
figure;
plot(dxVals, areas, '-s', 'LineWidth', 2, 'MarkerFaceColor', 'w')
set(gca, 'XDir', 'reverse')  % Invert x-axis
set(gca, 'XScale', 'log')    % Log scale
xlabel('Δx')
ylabel('Area (estimate of definite integral)')