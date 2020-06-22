function P = pwlinearFD(n,param)
% P = pwlinearFD(n,param)
% Compute the piecewise linear Fundamental Diagram function
%
% Nov 2019 - Guilhem Mariotte
%
% INPUTS
%---- n     : scalar or vector, accumulation [veh]
%---- param : Npts-by-2 matrix, points to define the piecewise linear FD
%             param(:,1) contains the ni of branch i start point [veh]
%             param(:,2) contains the Pi of branch i start point [veh.m/s]
%
% OUTPUTS
%---- P : scalar or vector (same size as n), production [veh.m/s]

ni = param(:,1); % abscissa of branch start point [veh]
Pi = param(:,2); % ordinate of branch start point [veh.m/s]

Npts = length(ni);
wi = zeros(Npts-1,1);
Pi0 = zeros(Npts-1,1);
for i = 1:(Npts-1)
    wi(i) = (Pi(i+1) - Pi(i)) / (ni(i+1) - ni(i));
    Pi0(i) = (Pi(i)*ni(i+1) - Pi(i+1)*ni(i)) / (ni(i+1) - ni(i));
end

wi = wi*ones(1,length(n));
Pi0 = Pi0*ones(1,length(n));
n = ones(size(wi,1),1)*n;

P = min(wi.*n + Pi0);
P = (0 <= P).*P;

end