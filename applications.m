%% Application of L'Hopital's:  Racing functions to infinitiy
% Which function grows faster to inf:
% to determine we can ratio them as
% lim x->inf f(x) / g(x) .

% Possible answers:
% inf = Numerator faster
% 0 = Denominator faster
% r = growth ratio settles to a value
% !!! Use L-Hopital's again and again to solve them!


% Note: The natural exponent grows to infinity faster than
% any power (or any non-infinitely differentiable function)


%% Second Derivative Test
% Purpose: distinguish minimum and maximum as critical points

% if ~x is a critical point of f, and
% if f'(~x) and f''(~x) exist, then:

% f''(~x) > 0 -> f(~x) is a local min
% f''(~x) < 0 -> f(~x) is a local max
% f''(~x) = 0 -> inconclusive

%% Compute and draw until second derivative
clear; clc; close all;

syms x
f = 2*pi*x + sin(pi*x) ; % we can test different functions here

secDerTest(f, x)

function secDerTest(f, x)
    df1 = diff(f, x);
    df2 = diff(f, x, 2);
    
    df1_zeros = double(solve(df1 == 0, x, 'Real', true));
    
    if isempty(df1_zeros)
        warning('No critical points found!')
        return
    end
    
    f_str = ['f = ' char(f)];
    df1_str = ['f'' = ' char(df1)];
    df2_str = ['f'''' = ' char(df2)];
    
    f_at_zeros = double(subs(f, x, df1_zeros));      % Eklendi
    df2_at_zeros = double(subs(df2, x, df1_zeros));
    
    figure(1);
    fplot(f, 'LineWidth', 2)
    hold on
    fplot(df1, 'LineWidth', 2)
    fplot(df2, 'LineWidth', 2)
    ylim([-40 40])
    grid on
    yline(0, '--k', 'LineWidth', 1)
    plot(df1_zeros, f_at_zeros, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
    plot(df1_zeros, df2_at_zeros, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', [0.5 0 0.5])
    hold off
    legend(f_str, df1_str, df2_str, 'y=0', 'f''=0', 'f''''=0')
    
    for i = 1:length(df1_zeros)
        xi = df1_zeros(i);
        df2_val = double(subs(df2, x, xi));  % double() ile çevir
        
        fprintf('Point %d: x = %.4f\n', i, xi);
        if df2_val > 0
            fprintf('  Result: Local MINIMUM\n\n');
        elseif df2_val < 0
            fprintf('  Result: Local MAXIMUM\n\n');
        else
            fprintf('  Result: INCONCLUSIVE\n\n');
        end
    end
end


%% Linear Approximations
% finding the value of the function x0 is the goal.
% find a point near the target: point a.
% compute the derivative at point a. extend that derivative to the x0.
% it is a better estimate compared to the near target a = x0...

% formula: L(x0) = f(a) + f'(a)*(x0-a)
% L is linear approximation
% intercept = f(a)
% slope = f'(a) 
% dx = x0-a (distance)

% This might fail due to the distance (dx). But it depends on the function.
% What is acceptable error is also depends on the field.
% Approximation over or under estimate? Depends on the curvature, we can do
%   second derivative test.


