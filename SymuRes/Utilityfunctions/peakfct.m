function peakf = peakfct(t,tstartpeak,tendpeak,peakfactor)
% peakf = peakfct(t,tstartpeak,tendpeak,peakfactor)
% Return a peak function (gaussian model)
%
% INPUTS
%---- t          : scalar or vector, time [s]
%---- tstartpeak : scalar, peak start time [s]
%---- tendpeak   : scalar, peak end time [s]
%---- peakfactor : scalar, increase factor (multiplier), amplitude of the peak
%
% OUTPUTS
%--- peakf : vector, same size as t, function values

tpeak = (tendpeak + tstartpeak)/2;
peakwidth = (tendpeak - tstartpeak)/8;

peakf = 1 + (tstartpeak <= t).*(t <= tendpeak).*(peakfactor - 1).*exp(-(t - tpeak).^2./(2*peakwidth^2));

end