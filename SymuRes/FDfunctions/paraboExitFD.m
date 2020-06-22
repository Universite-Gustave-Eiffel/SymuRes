function P = paraboExitFD(n,param)
% P = paraboExitFD(n,param)
% Compute the parabolic exit function (from the parabolic FD model)
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

P = (0 <= n).*(n <= nc).*paraboFD(n,param) + (nc < n).*(n < nj).*Pc;

end