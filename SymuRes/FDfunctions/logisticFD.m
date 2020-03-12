function q = logisticFD(k,param)
% q = logisticFD(k,param)
% Compute the 4 parameter logistic Fundamental Diagram function
% Wang et al. (2010)
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
a = (kj - kc)/8; % logistic range parameter
b = 1; % logistic shape parameter (lopsidedness)

q = (k < kj).*u.*k./((1 + exp((k - kc)./a)).^b);

end