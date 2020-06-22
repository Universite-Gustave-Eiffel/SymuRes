function P = semiparaboFD(n,param)
% P = semiparaboFD(n,param)
% Compute the semi-parabolic Fundamental Diagram function
% (n < nc: parabolic part, n > nc: linear part)
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

w = Pc/(nj - nc); % congestion wave speed [m/s]

P = (0 <= n).*(n <= nc).* (Pc/nc^2.*n.*(2*nc - n)) + ...
    (nc < n).*(n < nj).*  (Pc - w.*(n - nc));

end