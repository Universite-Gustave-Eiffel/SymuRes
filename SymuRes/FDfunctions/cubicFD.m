function P = cubicFD(n,param)
% P = cubicFD(n,param)
% Compute a third-order fundamental diagram
% Valid only for 1/3*nj < nc < 2/3*nj, and nc must be different from nj/2
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

b = nj;
c = (3*nc - 2*nj)./(2 - nj./nc);
a = Pc./(nc.*(nc - b).*(nc - c));

P = (0 <= n).*(n < nj).*a.*n.*(n - b).*(n - c);

end