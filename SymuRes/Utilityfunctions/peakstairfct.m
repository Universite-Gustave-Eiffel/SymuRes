function peakf = peakstairfct(t,tstartpeak,tendpeak,peakvalue,Nstairs)
% peakf = peakstairfct(t,tstartpeak,tendpeak,peakvalue,Nstairs)
% Return a peak function (gaussian model) by stairs
%
% INPUTS
%---- t          : scalar or vector, time [s]
%---- tstartpeak : scalar, peak start time [s]
%---- tendpeak   : scalar, peak end time [s]
%---- peakfactor : scalar, increase factor (multiplier), amplitude of the peak
%---- Nstairs    : scalar, number of stairs
%
% OUTPUTS
%--- peakf : vector, same size as t, stair function

Nt = length(t);
Npts = 1000; % number of points to calculate the mean value for a stair

timestair = linspace(tstartpeak,tendpeak,Nstairs+1); % time of each stair edges

peakf = zeros(1,Nt);

for i = 1:Nstairs
    stairt = linspace(timestair(i),timestair(i+1),Npts);
    stairvalue = mean(peakfct(stairt,tstartpeak,tendpeak,peakvalue));
    peakf = peakf + (timestair(i) <= t).*(t < timestair(i+1)).*stairvalue;
end

peakf = peakf + (t < tstartpeak) + (tendpeak <= t); % default value = 1 if t is not in the peak interval

end