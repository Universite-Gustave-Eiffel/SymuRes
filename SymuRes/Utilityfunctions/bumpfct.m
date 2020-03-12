function bumpf = bumpfct(t,tstartpeak,tendpeak,peakfactor1,peakfactor2)
% bumpf = bumpfct(t,tstartpeak,tendpeak,peakfactor1,peakfactor2)
% Return a bump function (sinusoidal model): a peak function with final
% value different from initial value
%
% INPUTS
%---- t           : scalar or vector, time [s]
%---- tstartpeak  : scalar, peak start time [s]
%---- tendpeak    : scalar, peak end time [s]
%---- peakfactor1 : scalar, multiplier factor, amplitude of the peak
%---- peakfactor2 : scalar, multiplier factor, amplitude of the new value
%                   after the peak
%
% OUTPUTS
%--- bumpf : vector, same size as t, function values

peakwidth = tendpeak - tstartpeak;
tmidpeak = (tstartpeak + tendpeak)/2;

bumpf = 1 + (tstartpeak <= t).*(t < tmidpeak).*(peakfactor1 - 1).*(1 + sin(2*pi*(t - tstartpeak)./peakwidth - pi/2))./2 ...
    + (tmidpeak <= t).*(t < tendpeak).*(peakfactor2 - 1 + (peakfactor1 - peakfactor2).*(1 + sin(2*pi*(t - tstartpeak)./peakwidth - pi/2))./2) ...
    + (tendpeak <= t).*(peakfactor2 - 1);

end