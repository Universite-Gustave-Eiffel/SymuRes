function P = expoFD(n,param)
% P = expoFD(n,param)
% Compute the exponential Fundamental Diagram function (Drake et al., 1967)
% Northwestern Model (1967)
% Only nc and Pc are required to define its shape
%
% Nov 2019 - Guilhem Mariotte
%
% INPUTS
%---- n     : scalar or vector, accumulation [veh], range of n must be between 0 and nj
%---- param : vector, parameters of the FD function = [nj nc Pc]
%
% OUTPUTS
%---- P : scalar or vector (same size as n), production [veh.m/s]

nj = param(1); % jam accumulation (max. accumulation) [veh]
nc = param(2); % critical accumulation [veh]
Pc = param(3); % critical production (max. production) [veh.m/s]

u = Pc/nc; % free-flow speed [m/s]

P = (0 <= n).*(n < nj).*u.*n.*exp(-(n./nc).^2/2);

end