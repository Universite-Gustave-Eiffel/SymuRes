function q = trapzFD(k,param)
% q = trapzFD(k,param)
% Compute the trapezoid-shaped Fundamental Diagram function 
%
% INPUTS
%---- k     : scalar or vector, vehicle density [veh/m], range of k must be between 0 and kj
%---- param : vector, parameters of the FD function = [kj kc1 kc2 qc1 qc2]
%
% OUTPUTS
%---- q : scalar or vector (same size as k), vehicle flow [veh/s]

kj = param(1); % jam density (max. density) [veh/m]
kc1 = param(2); % first critical density [veh/m]
kc2 = param(3); % second critical density [veh/m]
qc1 = param(4); % first critical flow [veh/s]
qc2 = param(5); % second critical flow [veh/s]

u = qc1/kc1; % free-flow speed [m/s]
w = qc2/(kj - kc2); % congestion wave speed [m/s]

q = (k <= kc1).*(u.*k) + (k > kc1).*(k <= kc2).*(qc1 + (qc2 - qc1)/(kc2 - kc1).*(k - kc1)) + (k > kc2).*(qc2 - w.*(k - kc2));

end