function [tval, fval] = fstair(f,tstart,tend,Nstairs)
% [tval, fval] = fstair(f,tstart,tend,Nstairs)
% Return function values by stairs
%
% INPUTS
%---- f       : function handle
%---- tstart  : scalar, start time [s]
%---- tend    : scalar, end time [s]
%---- Nstairs : scalar, number of stairs
%
% OUTPUTS
%---- tval : vector, size Nstairs+1, time values at the stair edges
%---- fval : vector, size Nstairs+1, stair function values

stairwidth = (tend - tstart)/Nstairs; % width of each stair

Npts = 1000; % number of points to calculate the mean value for a stair

tval = zeros(1,Nstairs+1);
fval = zeros(1,Nstairs+1);

tval(1) = tstart;
for i = 1:Nstairs
    tval(i+1) = tval(i) + stairwidth;
    
    stairt = linspace(tval(i),tval(i+1),Npts);
    stairvalue = mean(f(stairt));
    fval(i) = stairvalue;
end

fval(Nstairs+1) = stairvalue;

end