function P = logisticFD(n,param)
% P = logisticFD(n,param)
% Compute the 4 parameter logistic Fundamental Diagram function
% Wang et al. (2010)
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
a = (nj - nc)/8; % logistic range parameter
b = 1; % logistic shape parameter (lopsidedness)

P = (0 <= n).*(n < nj).*u.*n./((1 + exp((n - nc)./a)).^b);

end