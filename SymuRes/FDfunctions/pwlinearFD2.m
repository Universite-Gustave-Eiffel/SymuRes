function q = pwlinearFD2(k,param)
% q = pwlinearFD2(k,param)
% Compute the piecewise linear Fundamental Diagram function
%
% INPUTS
%---- k     : scalar or vector, vehicle density [veh/m], range of k must be between 0 and kj
%---- param : Npts-by-2 matrix, points to define the piecewise linear FD
%             column 1 contains the slopes [m/s]
%             column 2 contains the k-intercepts corresponding to each slope [veh/m]
%
% OUTPUTS
%---- q : scalar or vector (same size as k), vehicle flow [veh/s]

wc = param(:,1); % slope of the FD branches [m/s]
kjc = param(:,2); % k-intercept of the FD branches [veh/m]

wc = wc*ones(1,length(k));
kjc = kjc*ones(1,length(k));
k = ones(length(wc),1)*k;

q = min(wc.*(k - kjc));