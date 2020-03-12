function q = expoFD(k,param)
% q = expoFD(k,param)
% Compute the exponential Fundamental Diagram function (Drake et al., 1967)
% Northwestern Model (1967)
% !! in fact only kc and qc are required to define its shape
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

u = qc/kc; % free-flow speed [m/s]

q = (k < kj).*u.*k.*exp(-(k./kc).^2/2);

end