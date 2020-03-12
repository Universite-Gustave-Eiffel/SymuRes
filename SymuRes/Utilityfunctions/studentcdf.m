function p = studentcdf(x0,ddl,numpts,xmin)
% p = studentcdf(x0,ddl,numpts,xmin)
% Compute the Student's cumulative distribution function
% (equivalent to the Matlab function 'tcdf' from the Statistics Toolbox)
%
% INPUTS
%---- x0     : scalar or vector, value for which the probability is calculated
%---- ddl    : integer, degree of freedom of the Student distribution
%---- numpts : integer, number of points for the calculation of the discrete integral
%---- xmin   : scalar, min. bound of the integral (theoretically equal to -Inf)
%
% OUTPUTS
%---- p : scalar or vector, probability given by the Student's CDF at the value x0

x = zeros(numpts,length(x0));
for i = 1:numpts
    x(i,:) = xmin + (x0 - xmin).*(i - 1)./(numpts - 1);
end
dx = (x0 - xmin)./(numpts - 1);

if ddl < 300 % numerical errors for high values of degree of freedom
    p = dx.*gamma((ddl+1)/2)/gamma(ddl/2)*1/sqrt(pi*ddl).*sum(1./(1 + x.^2./ddl).^((ddl+1)/2));
else % approximation by a normal distribution (mu = 0, sigma = 1)
    p = dx.*1./sqrt(2*pi).*sum(exp(-x.^2./2));
end

end