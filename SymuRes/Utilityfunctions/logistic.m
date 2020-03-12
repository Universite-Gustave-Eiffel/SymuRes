function y = logistic(x,K,a,lambda)
% y = logistic(x,K,a,lambda)
% Compute the logistic function
%
% INPUTS
%---- x      : scalar or vector, points where evaluating the function
%---- K      : scalar, scale parameter (y-range)
%---- a      : scalar > 0, x-range
%---- lambda : scalar > 0, delay parameter
%
% OUTPUTS
%---- y : values of the function

r = 2*6/a;
x0 = 6/r;

y = K./(1 + lambda.*exp(-r.*(x - x0)));

end