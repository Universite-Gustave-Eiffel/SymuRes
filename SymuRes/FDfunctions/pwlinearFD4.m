function q = pwlinearFD4(k,param)
% q = pwlinearFD4(k,param)
% Compute the piecewise linear Fundamental Diagram function, with inputs as
% slope-intercept values (wc,qc0).
%
% INPUTS
%---- k     : scalar or vector, vehicle density [veh/m], range of k must be between 0 and kj
%---- param : Npts-by-2 matrix, points to define the piecewise linear FD
%             column 1 contains the slopes wc [m/s]
%             column 2 contains the q-intercepts qc0 corresponding to each slope [veh/s]
%
% OUTPUTS
%---- q : scalar or vector (same size as k), vehicle flow [veh/s]

wc = param(:,1); % slope of the FD branches [m/s]
qc0 = param(:,2); % q-intercept of the FD branches [veh/s]

wc = wc*ones(1,length(k));
qc0 = qc0*ones(1,length(k));
k = ones(size(wc,1),1)*k;

q = min(wc.*k + qc0);

end