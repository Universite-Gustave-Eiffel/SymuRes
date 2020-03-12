function [tval, fval] = fstair2(f,tstart,tend,Nstairs)
% [tval, fval] = fstair2(f,tstart,tend,Nstairs)
% Return function values by stairs
%
% INPUTS
%---- f       : function handle
%---- tstart  : scalar, start time [s]
%---- tend    : scalar, end time [s]
%---- Nstairs : scalar, number of stairs
%
% OUTPUTS
%--- tval : vector, size 2*Nstairs, time values at the stair edges
%--- fval : vector, size 2*Nstairs, stair function values

stairwidth = (tend - tstart)/Nstairs; % width of each stair

margin = 0.001; % margin between two time values at a stair interface
Npts = 1000; % number of points to calculate the mean value for a stair

tval = zeros(1,2*Nstairs);
fval = zeros(1,2*Nstairs);

i = 1;
tval(1+2*(i-1)) = tstart;
tval(2*i) = (1 - margin)*(tstart + stairwidth);

stairt = linspace(tval(1+2*(i-1)),tval(2*i),Npts);
stairvalue = mean(f(stairt));
fval(1+2*(i-1)) = stairvalue;
fval(2*i) = stairvalue;

for i = 2:Nstairs
    t0 = tval(2*(i-1)-1) + stairwidth;
    tval(1+2*(i-1)) = t0;
    tval(2*i) = (1 - margin)*(t0 + stairwidth);
    
    stairt = linspace(tval(1+2*(i-1)),tval(2*i),Npts);
    stairvalue = mean(f(stairt));
    fval(1+2*(i-1)) = stairvalue;
    fval(2*i) = stairvalue;
end

end