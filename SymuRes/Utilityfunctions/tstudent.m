function t0 = tstudent(p0,ddl)
% t0 = tstudent(p0,ddl)
% Compute the inverse of the Student's cumulative distribution function
% (equivalent to the Matlab function 'tinv' from the Statistics Toolbox)
%
% INPUTS
%---- p0  : scalar or vector, probability in [0,1] given by the Student's CDF at t0
%---- ddl : degree of freedom of the Student distribution
%
% OUTPUTS
%---- t0 : scalar or vector, Student's t value which gives the probability p0

Npts = length(p0);

t0 = zeros(Npts);

for i = 1:Npts
    f = @(x) studentcdf(x,ddl,1e5,-1e3);
    t0(i) = findintersec(f,p0(i),0,50,1e-5,0.01);
end

end