function [t, f] = pwfct(tval,fval,Dt,tstart,tend)
% [t, f] = pwfct(tval,fval,Dt,tstart,tend)
% Return the values of a piecewise linear function evaluated at regular
% sampled points with a given time step
%
% INPUTS
%---- tval   : vector, time values that define the function points
%---- fval   : vector, same size as tval, function values
%---- Dt     : scalar, time step
%---- tstart : scalar, start time
%---- tend   : scalar, end time
%
% OUTPUTS
%---- t : vector, sampled times with time step Dt
%---- f : vector, same size as t, stair function values

Nt = floor((tend - tstart)/Dt) + 1; % number of times
t = tstart:Dt:tend; % vector of times

f = zeros(1,Nt);

Nval = length(fval); % number of stairs

i = 1;
f = f + (t < tval(i)).*fval(i);

for i = 2:Nval
    f = f + (tval(i-1) <= t).*(t < tval(i)).*(fval(i-1) + (fval(i) - fval(i-1))./(tval(i) - tval(i-1)).*(t - tval(i-1)));
end

i = Nval;
f = f + (tval(i) <= t).*fval(i);

end