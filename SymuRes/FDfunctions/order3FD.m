function q = order3FD(k,param)
% q = order3FD(k,param)
% Compute a third order Fundamental Diagram function
% !! in fact only kj and qc are required to define its shape (kc = kj/3)
%
% INPUTS
%---- k     : scalar or vector, vehicle density [veh/m], range of k must be between 0 and kj
%---- param : vector, parameters of the FD function = [kj kc qc]
%
% OUTPUTS
%---- q : scalar or vector (same size as k), vehicle flow [veh/s]

kj = param(1); % jam density (max. density) [veh/m]
kc = kj/3; % critical density [veh/m]
qc = param(2); % critical flow (max. flow) [veh/s]

u = qc/kc; % free-flow speed [m/s]

% q = 9/4*k.*u.*(1 - k./kj).^2;
q = k.*u.*(1 - k./kj).^2;

end