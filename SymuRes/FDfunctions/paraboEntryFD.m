function P = paraboEntryFD(n,param)
% P = paraboFD(n,param)
% Compute the parabolic entry function (from the parabolic FD model)
%
% Nov 2019 - Guilhem Mariotte
%
% INPUTS
%---- n     : scalar or vector, accumulation [veh], range of n must be between 0 and nj
%---- param : vector, parameters of the entry function = [nj nc Pc nc1 nc2 Pc1]
%
% OUTPUTS
%---- P : scalar or vector (same size as n), production [veh.m/s]

FDparam = param(1:3); % FD function parameters = [nj nc Pc]
nj = FDparam(1);

nc1 = param(4); % critical accumulation for free-flow: 0 < nc1 <= nc [veh]
nc2 = param(5); % critical accumulation for the transition to congestion: nc1 <= nc <= nc2 [veh]
Pc1 = param(6); % free-flow max production: Pc1 >= Pc [veh.m/s]
Pc2 = paraboFD(nc2,FDparam);

P = (n <= nc1).*            Pc1 + ...
    (nc1 < n).*(n <= nc2).* (Pc1 + (n - nc1)./(nc2 - nc1).*(Pc2 - Pc1)) + ...
    (nc2 < n).*(n < nj).*   paraboFD(n,FDparam);

end