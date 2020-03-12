function peakf = peakfct2(t,tstartpeak,tendpeak,peakfactor)
% peakf = peakfct2(t,tstartpeak,tendpeak,peakfactor)
% Return a peak function (sinusoidal model)
%
% INPUTS
%---- t          : scalar or vector, time [s]
%---- tstartpeak : scalar, peak start time [s]
%---- tendpeak   : scalar, peak end time [s]
%---- peakfactor : scalar, increase factor (multiplier), amplitude of the peak
%
% OUTPUTS
%--- peakf : vector, same size as t, function values

peakwidth = tendpeak - tstartpeak;

peakf = 1 + (tstartpeak <= t).*(t <= tendpeak).*(peakfactor - 1).*(1 + sin(2*pi*(t - tstartpeak)./peakwidth - pi/2))./2;

end