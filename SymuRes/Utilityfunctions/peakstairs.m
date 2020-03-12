function [peakt, peakf] = peakstairs(t,tstartpeak,tendpeak,peakvalue,Nstairs)
% [peakt, peakf] = peakstairs(t,tstartpeak,tendpeak,peakvalue,Nstairs)
% Return a peak function (gaussian model) by stairs
%
% INPUTS
%---- t          : vector, time [s]
%---- tstartpeak : scalar, peak start time [s]
%---- tendpeak   : scalar, peak end time [s]
%---- peakfactor : scalar, increase factor (multiplier), amplitude of the peak
%---- Nstairs    : scalar, number of stairs
%
% OUTPUTS
%--- peakt : vector, size 2*(Nstairs + 2), time values at the stair edges
%--- peakf : vector, size 2*(Nstairs + 2), stair function values

stairwidth = (tendpeak - tstartpeak)/Nstairs; % width of each stair

margin = 0.001; % margin between two time values at a stair interface
Npts = 1000; % number of points to calculate the mean value for a stair

peakt = zeros(1,2*Nstairs);
peakf = zeros(1,2*Nstairs);

i = 1;
peakt(1+2*(i-1)) = tstartpeak;
peakt(2*i) = (1 - margin)*(tstartpeak + stairwidth);

stairt = linspace(peakt(1+2*(i-1)),peakt(2*i),Npts);
stairvalue = mean(peakfct(stairt,tstartpeak,tendpeak,peakvalue));
peakf(1+2*(i-1)) = stairvalue;
peakf(2*i) = stairvalue;

for i = 2:Nstairs
    tvalue = peakt(2*(i-1)-1) + stairwidth;
    peakt(1+2*(i-1)) = tvalue;
    peakt(2*i) = (1 - margin)*(tvalue + stairwidth);
    
    stairt = linspace(peakt(1+2*(i-1)),peakt(2*i),Npts);
    stairvalue = mean(peakfct(stairt,tstartpeak,tendpeak,peakvalue));
    peakf(1+2*(i-1)) = stairvalue;
    peakf(2*i) = stairvalue;
end

% Add the initial and final time points
Nt = length(t);
peakt = [t(1) (1 - margin)*tstartpeak peakt tendpeak t(Nt)];
peakf = [1 1 peakf 1 1];

end