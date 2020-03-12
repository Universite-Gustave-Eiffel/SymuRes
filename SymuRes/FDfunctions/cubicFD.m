function q = cubicFD(k,param)
% q = cubicFD(k,param)
% Compute a third-order fundamental diagram
% Valid only for 1/3*kj < kc < 2/3*kj, and kc must be different from kj/2
%
% INPUTS
%---- k     : scalar or vector, vehicle density [veh/m], range of k must be between 0 and kj
%---- param : vector, parameters of the FD function = [kj kc qc]
%
% OUTPUTS
%---- q : scalar or vector (same size as k), vehicle flow [veh/s]

kj = param(1); % jam density (max. density) [veh/m]
kc = param(2); % critical density [veh/m]
qc = param(3); % critical flow (max. flow) [veh/s]

b = kj;
c = (3*kc - 2*kj)./(2 - kj./kc);
a = qc./(kc.*(kc - b).*(kc - c));

q = a.*k.*(k - b).*(k - c);

end