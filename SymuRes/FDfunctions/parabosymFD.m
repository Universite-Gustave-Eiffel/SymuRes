function P = parabosymFD(n,param)
% P = parabosymFD(n,param)
% Compute the parabolic Fundamental Diagram function (symmetrical)
% Greenshields Model (1935)
% Only nj and Pc are required to define its shape (nc = nj/2)
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
nc = nj/2; % critical accumulation [veh]
Pc = param(3); % critical production (max. production) [veh.m/s]

u = 2*Pc/nc; % free-flow speed [m/s]

P = (0 <= n).*(n < nj).*u.*n.*(1 - n./nj);

end