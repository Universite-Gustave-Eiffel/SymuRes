function trayf = trayfct(t,tstart,ttray1,ttray2,tend,trayfactor,finalfactor)
% trayf = trayfct(t,tstart,ttray1,ttray2,tend,trayfactor,finalfactor)
% Return a tray function: initial value, then sinusoidal transition to a
% steady value, and finally another sinusoidal transition to a final value
%
% INPUTS
%---- t           : scalar or vector, time [s]
%---- tstart      : scalar, tray start time [s]
%---- ttray1      : scalar, tray start transition end time [s]
%---- ttray2      : scalar, tray end transition start time [s]
%---- tend        : scalar, tray end time [s]
%---- trayfactor  : scalar, multiplier factor, amplitude of the tray
%---- finalfactor : scalar, multiplier factor, amplitude of the new value
%                   after the tray
%
% OUTPUTS
%--- bumpf : vector, same size as t, function values

traywidth1 = 2*(ttray1 - tstart);
traywidth2 = 2*(tend - ttray2);

trayf = 1 + (tstart <= t).*(t < ttray1).*(trayfactor - 1).*(1 + sin(2*pi*(t - tstart)./traywidth1 - pi/2))./2 ...
    + (ttray1 <= t).*(t < ttray2).*(trayfactor - 1) ...
    + (ttray2 <= t).*(t < tend).*(finalfactor - 1 + (trayfactor - finalfactor).*(1 + sin(2*pi*(t - tend + traywidth2)./traywidth2 - pi/2))./2) ...
    + (tend <= t).*(finalfactor - 1);

end