%% Linearity of differentiation
% Why it is linear? Proof:
% k(x) = a*f(x) + b*g(x) 

% k' = lim h -> 0 ( k(x+h) - k(x) ) / h
%    = lim h -> 0 ( (a*f(x+h) + b*g(x+h)) - a*f(x) - b*g(x) ) / h
%    = lim h -> 0 ( a*(f(x+h) - f(x)) + b*(g(x+h) - g(x)) ) / h
% k' = ( a* f' ) + ( b* g')

%% THEOREM:  Differentiability implies continuity 
% If derivative exist at point x0,
% then the function is continuous at x0.

% We start with definition of continuity: lim x->a f(x) = f(a)

% lim x->a (f(x) - f(a)) = 0  % we can do it because there is no x in f(a)
% its a constant.

% Then we multiply and divide by x-a..
% lim x->a (f(x) - f(a)) / (x - a) * lim x->a (x - a) = 0
% First part is derivative of f and second part is x = a since x -> a and 0
% f'(x) * 0 = 0 : CORRECT !

% OPPOSITE WAY IS NOT CORRECT! IF FUNCTION IS CONTINOUS IT DOES NOT MEAN IT
% HAS TO HAVE DERIVATIVE like f(x) = |x|.

%% PRODUCT RULE
% ( f(x)*g(x) )' = f'(x)*g(x) + f(x)*g'(x) 








