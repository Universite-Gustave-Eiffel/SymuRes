function xres = findintersec(f,yintersec,x0,maxiter,crit,eps)
% xres = findintersec(f,yintersec,x0,maxiter,crit,eps)
% Find the intersection xres of the function f(x) with a constant value y
% (Equivalent to find the zero xres of f(x)-y using Newton's method)
%
% INPUTS
%---- f         : function for which the intersection has to be found
%---- yintersec : scalar, y value of the intersection
%---- x0        : scalar, start point of the algorithm
%---- maxiter   : int, maximum number of iterations
%---- crit      : scalar, convergence criterium, minimum acceptable
%                 difference between two successive iteration
%---- eps       : scalar, espilon to calculate the derivative of f(x)
%
% OUTPUTS
%---- xres : scalar, result of the algorithm, f(xres) = yintersec

f0 = @(x) f(x) - yintersec;

xres = x0;
xres1 = x0 + 2*crit;
niter = 1;
while abs(xres - xres1) > crit && niter < maxiter
    xres1 = xres;
    df = (f0(xres + eps) - f0(xres - eps))/(2*eps);
    xres = xres - f0(xres)/df;
    niter = niter + 1;
end

end