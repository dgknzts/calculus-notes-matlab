%% Probability Density Functions (PDFs)
% PDF is a function that maps events onto probabilities.
% Two constrains:
% 1) p(x) >= 0, All x in X  (Non-Negativity)
% 2) integral -inf to inf p(x) dx = 1 (Normalization)

%% Normal (Gaussian) Distribution
syms mu sigma x

% PDF formula
pdf_normal = (1/(sigma*sqrt(2*pi))) * exp(-(x-mu)^2 / (2*sigma^2));

% Substitute values
mu_val = 0;
sigma_val = 2;
pdf_normal_vals = subs(pdf_normal, [mu, sigma], [mu_val, sigma_val]);

figure(1);
fplot(pdf_normal_vals, [-5, 5], 'LineWidth', 2);
xlabel('x'); ylabel('f(x)');
title(sprintf('Normal Distribution (μ=%d, σ=%d)', mu_val, sigma_val));
grid on;

% Convert symbolic to function handle
pdf_func = matlabFunction(pdf_normal_vals);

% Now evaluate over the range
xx = linspace(-20, 20, 500); 
yy = pdf_func(xx);

% Sum
sum_val = sum(yy) * (xx(2) - xx(1));
fprintf('Sum: %f\n', sum_val);

%% Logistic Distribution PDF
clear; clc; close all;

syms x mu s
assume(s > 0);

% Logistic PDF formula
pdf_logistic = exp(-(x-mu)/s) / (s * (1 + exp(-(x-mu)/s))^2);


% Substitute parameter values
mu_val = 0;   % location (mean)
s_val = 1;    % scale (controls spread)

pdf_logistic_vals = subs(pdf_logistic, [mu, s], [mu_val, s_val]);

% Verify it integrates to 1
total_prob = int(pdf_logistic_vals, x, -inf, inf);
fprintf('\nTotal probability: %s (should be 1)\n', char(total_prob));

% Plot
figure(1);
fplot(pdf_logistic_vals, [-10, 10], 'LineWidth', 2);
xlabel('x'); ylabel('f(x)');
title(sprintf('Logistic PDF (μ=%d, s=%d)', mu_val, s_val));
grid on;

%% Exercise 3: Exact definite integral
clear; clc; close all;

syms x

% Define the integrand (logistic PDF with specific parameters)
pdf = exp(-x/2) / (2*(1 + exp(-x/2))^2);

% Compute the exact definite integral from 0 to π
area_exact = int(pdf, x, 0, pi);

fprintf('Exact integral:\n');
pretty(area_exact);
fprintf('\nNumerical value: %.6f\n', double(area_exact));

% Create numerical version for plotting
pdf_func = matlabFunction(pdf);
xx = linspace(-5, 20, 1000);
yy = pdf_func(xx);

% Find the shaded region (0 to π)
x_region = linspace(0, pi, 100);
y_region = pdf_func(x_region);

% Plot
figure(1);
plot(xx, yy, 'r', 'LineWidth', 2);
hold on;

% Fill the area from 0 to π
fill([x_region, fliplr(x_region)], [y_region, zeros(size(y_region))], ...
    'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Add vertical lines at boundaries
yline(0, 'k', 'LineWidth', 0.5);
plot([0, 0], [0, pdf_func(0)], 'k--', 'LineWidth', 1);
plot([pi, pi], [0, pdf_func(pi)], 'k--', 'LineWidth', 1);

xlabel('x');
ylabel('pdf');
title('Exact definite integral');
legend('Logistic pdf', sprintf('Area = %.3f', double(area_exact)));
grid on;
xlim([-5, 20]);
hold off;

% Display the formula
fprintf('\nIntegral formula:\n');
fprintf('∫₀^π [e^(-x/2) / (2(1 + e^(-x/2))²)] dx = -1/2 + 1/(e^(π/2) + 1)\n');

%% Exercise: Logistic CDF (Three Ways)
clear; clc; close all;

N = 31;
x = linspace(-5, 5, N);
dx = x(2) - x(1);

% Logistic PDF formula
pdf_values = exp(-x) ./ (1 + exp(-x)).^2;

% Method 1: Analytical CDF
cdf_analytical = 1 ./ (1 + exp(-x));

% Method 2: Cumulative sum
cdf_cumsum = cumsum(pdf_values) * dx;

% Method 3: Numerical integration
cdf_integral = zeros(1, N);
for i = 1:N
    cdf_integral(i) = integral(@(t) exp(-t)./(1+exp(-t)).^2, -5, x(i));
end

% Plot
figure(1);
subplot(2, 1, 1);
plot(x, pdf_values, 'b-', 'LineWidth', 2);
xlabel('x'); ylabel('pdf');
title(sprintf('Logistic PDF (N = %d)', N));
grid on;

subplot(2, 1, 2);
hold on;
plot(x, cdf_analytical, 'b-',  'LineWidth', 2, 'DisplayName', 'Analytical');
plot(x, cdf_cumsum,     'r--', 'LineWidth', 2, 'DisplayName', 'Cumulative sum');
plot(x, cdf_integral,   'g-',  'LineWidth', 2, 'DisplayName', 'Numerical integration');
hold off;
xlabel('x'); ylabel('cdf');
title(sprintf('Logistic CDF (N = %d)', N));
legend('Location', 'northwest');
grid on;

%% Exercise Area from PDF and CDF
clear; clc; close all;

mu = sqrt(2)/2;
sigma = 1/pi;
N = 301;
x = linspace(-2, 2, N);

pdf_values = normpdf(x, mu, sigma);
cdf_values = normcdf(x, mu, sigma);

% Shaded region: from -inf to sqrt(2)/pi
b = sqrt(2)/pi;

% Area from PDF using integral
area_from_pdf = integral(@(t) normpdf(t, mu, sigma), -inf, b);

% Area from CDF by indexing
[~, idx_b] = min(abs(x - b));
area_from_cdf = cdf_values(idx_b);

fprintf('Area from pdf: %.4f\n', area_from_pdf);
fprintf('Area from cdf: %.4f\n', area_from_cdf);

% Exact result using Symbolic Math Toolbox
syms u
c(u) = erf(sym(pi) * (u - sqrt(sym(2))/2) / sqrt(sym(2))) / 2 + sym(1)/2;
result_exact = double(c(sqrt(sym(2))/sym(pi)));
fprintf('Exact:         %.4f\n', result_exact);

% Plot
figure('Position', [100, 100, 600, 700]);

% PDF plot
subplot(2, 1, 1);
plot(x, pdf_values, 'm-', 'LineWidth', 2);
hold on;
% Shade the area from start to b
idx_shade = x <= b;
area(x(idx_shade), pdf_values(idx_shade), 'FaceColor', [0.6 0.6 1], ...
     'FaceAlpha', 0.4, 'EdgeColor', 'none');
hold off;
xlabel('x'); ylabel('Probability density');
title(sprintf('Normal pdf (\\mu = \\surd2/2, \\sigma = \\pi^{-1}, N = %d)', N));

% CDF plot
subplot(2, 1, 2);
plot(x, cdf_values, 'm-', 'LineWidth', 2);
hold on;
plot([b b], [0 cdf_values(idx_b)], 'b--', 'LineWidth', 1);
plot([-2 b], [cdf_values(idx_b) cdf_values(idx_b)], 'b--', 'LineWidth', 1);
hold off;
xlabel('x'); ylabel('Cumulative probability');
title('Normal cdf');