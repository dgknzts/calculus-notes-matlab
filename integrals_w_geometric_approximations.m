%% Riemann Sum
clear; clc; close all;
a = -0.5;
b = pi;
x = linspace(a, b, 15);
dx = x(2) - x(1); 
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

%% Net and Total Are
% Net area = Positive - negative values
% Total = Sum areas ignoring sign
% To calculate by hand: find roots -> calc areas -> sum parts

%% Trapezidal Rule
% Trapezoid Area = width * mean(h1,h2)
clear; clc; close all;
a = -0.6;
b = 1;
x = linspace(a, b, 7);
dx = x(2) - x(1); 
fx = x.^2 - 0.5;
[~, fx0idx] = min(abs(fx-0));

figure(1);
hold on
totArea = 0;
negArea = 0;
posArea = 0;
% Draw Riemann bars (left rule)
for i = 1:length(x)-1
    totArea = dx * ((fx(i) + fx(i+1)) / 2) + totArea;
    if i < fx0idx
    negArea = dx * ((fx(i) + fx(i+1)) / 2) + negArea;
    plot(polyshape([x(i), ...
        x(i),...
        x(i+1),...
        x(i+1)],...
        [0,fx(i),...
        fx(i+1),...
        0]),'FaceColor','red');
    else
    posArea = dx * ((fx(i) + fx(i+1)) / 2) + posArea;
    plot(polyshape([x(i), ...
        x(i),...
        x(i+1),...
        x(i+1)],...
        [0,fx(i),...
        fx(i+1),...
        0]),'FaceColor','green');
    end
end
netArea = posArea - negArea;
syms xa
fxa = xa^2 - 0.5;
fplot(fxa, [-0.6 1], 'k', 'LineWidth', 2)
xlim([-0.7, 1.1]);
xline(a, '--g', 'LineWidth', 1)
xline(b, '--r', 'LineWidth', 1)
yline(0, '--')
hold off
xlabel('x')
ylabel('y = f(x)')
title([sprintf('Total Area = %.2f \n', totArea), ...
    sprintf('Net Area = %.2f \n', netArea),...
    sprintf('deltaX = %.2f', dx)])

% Analytical Area
ifxa_negative = double(int(fxa, xa, a, x(fx0idx)));
ifxa_positive = double(int(fxa, xa, x(fx0idx), b));
net_analytical = ifxa_positive - ifxa_negative;
total_analytical = double(int(fxa, xa, a, b));
fprintf('Total Area      = %.2f\n', totArea)
fprintf('Total Area Anal = %.2f\n', total_analytical)
fprintf('Net Area        = %.2f\n', netArea)
fprintf('Net Area Anal   = %.2f\n', net_analytical)


%% Lebesgue Integrals
clear; clc; close all;
% Riemann : Partition the domain (x axis)
% Lebesgue : Partition the range (y axis) 
% Horizontal bars instead of vertical in lebek
% Compare Riemann vs Lebesgue

f = @(x) x.^2 + cos(2*pi*x) / 5;
a = 0; b = 1;

Ns = 4:40;  % Different partition numbers
riemann_results = zeros(size(Ns));
lebesgue_results = zeros(size(Ns));

for k = 1:length(Ns)
    n = Ns(k);
    
    % Riemann (partition x-axis)
    dx = (b-a) / n;
    riemann = 0;
    for i = 0:n-1
        x_mid = a + dx*i + dx/2;  % midpoint
        riemann = riemann + f(x_mid) * dx;  % height × width
    end
    riemann_results(k) = riemann;
    
    % Lebesgue (partition y-axis)
    x_fine = linspace(a, b, 1000);
    dx_fine = (b-a) / 1000;
    y_vals = f(x_fine);
    
    y_parts = linspace(min(y_vals), max(y_vals), n+1);  % y-axis partitions
    
    lebesgue = 0;
    for i = 1:n
        % Find x points where f(x) is in this y-band
        in_band = (y_vals >= y_parts(i)) & (y_vals < y_parts(i+1));
        
        measure = sum(in_band) * dx_fine;  % "width" of x values in this band
        avg_y = (y_parts(i) + y_parts(i+1)) / 2;  % average y in band
        
        lebesgue = lebesgue + avg_y * measure;
    end
    lebesgue_results(k) = lebesgue;
end

%% Plot
analytical = integral(f, a, b);

plot(Ns, riemann_results, 'b', 'LineWidth', 2)
hold on
plot(Ns, lebesgue_results, 'r', 'LineWidth', 2)
yline(analytical, '--k', 'LineWidth', 2)
hold off

xlabel('N (partitions)')
ylabel('Integral')
legend('Riemann', 'Lebesgue', 'Analytical')