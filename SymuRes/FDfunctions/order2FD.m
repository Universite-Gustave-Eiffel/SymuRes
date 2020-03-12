function q = order2FD(k,param)
% q = order2FD(k,param)
% Compute a second order Fundamental Diagram function
% !! in fact only kj and qc are required to define its shape (kc = kj/2)
%
% INPUTS
%---- k     : scalar or vector, vehicle density [veh/m], range of k must be between 0 and kj
%---- param : vector, parameters of the FD function = [kj kc qc]
%
% OUTPUTS
%---- q : scalar or vector (same size as k), vehicle flow [veh/s]

kj = param(1); % jam density (max. density) [veh/m]
kc = kj/2; % critical density [veh/m]
qc = param(2); % critical flow (max. flow) [veh/s]

u = 2*qc/kc; % free-flow speed [m/s]

q = k.*u.*(1 - k./kj);

end