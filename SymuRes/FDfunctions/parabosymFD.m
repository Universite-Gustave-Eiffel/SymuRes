function q = parabosymFD(k,param)
% q = parabosymFD(k,param)
% Compute the parabolic Fundamental Diagram function (symmetrical)
% Greenshields Model (1935)
% !! in fact only kj and qc are required to define its shape (kc = kj/2)
%
% INPUTS
%---- k     : scalar or vector, vehicle density [veh/m], range of k must be between 0 and kj
%---- param : vector, parameters of the FD function = [kj kc qc]
%
% OUTPUTS
%---- q : scalar or vector (same size as k), vehicle flow [veh/s]

kj = param(1); % jam density (max. density) [veh/m]
qc = param(3); % critical flow (max. flow) [veh/s]

q = 4*qc/kj^2*k.*(kj - k).^2;

end