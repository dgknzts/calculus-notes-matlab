%% Partial Fractions Algorithm
clear; clc;
syms x A B

fx = (5*x + 3) / (2*x^2 - 4*x - 6);
[num, den] = numden(fx);

% Get leading coefficient (coefficient of highest power)
coeffs_den = coeffs(den, x, 'All');
leading_coeff = coeffs_den(1); 

% Find roots of denominator
roots_den = solve(den == 0, x);
den_factor1 = leading_coeff * (x - roots_den(1));  % 2*(x+1)
den_factor2 = (x - roots_den(2));                   % (x-3)


f1 = A / den_factor1;
f2 = B / den_factor2;

% Get numerator when combined
[num_AB, ~] = numden(f1 + f2);

% Substitute x = roots_den(2) to eliminate A
B_ans = solve(subs(num, x, roots_den(2)) == subs(num_AB, x, roots_den(2)), B);

% Substitute x = roots_den(1) to eliminate B
A_ans = solve(subs(num, x, roots_den(1)) == subs(num_AB, x, roots_den(1)), A);

fprintf('A = %s\n', char(A_ans));
fprintf('B = %s\n', char(B_ans));

% Substitute A and B back into partial fractions
f1_final = subs(f1, A, A_ans);
f2_final = subs(f2, B, B_ans);

fprintf('\nPartial fractions:\n');
pretty(f1_final);
pretty(f2_final);

% Integrate each term separately
integral_f1 = int(f1_final, x);
integral_f2 = int(f2_final, x);

fprintf('\nTotal integral:\n');
total_integral = integral_f1 + integral_f2;
pretty(total_integral);