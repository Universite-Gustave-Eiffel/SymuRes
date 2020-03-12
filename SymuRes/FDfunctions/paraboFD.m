function q = paraboFD(k,param)
% q = paraboFD(k,param)
% Compute the parabolic Fundamental Diagram function (two parabolic arcs)
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

q = (k <= kc).*(qc./kc.^2.*k.*(2*kc - k)) + (k > kc).*(qc./(kj - kc).^2.*(kj - k).*(kj + k - 2*kc));

end