function q = pwlinearFD3(k,param)
% q = pwlinearFD3(k,param)
% Compute the piecewise linear Fundamental Diagram function with a scaling
% factor.
%
% INPUTS
%---- k     : scalar or vector, vehicle density [veh/m]
%---- param : 2-size cell
%             param{1} : scalar, apply a scaling factor in the function
%             param{2} : Npts-by-2 matrix, points to define the piecewise linear FD
%                 column 1 contains the kc of branch start point [veh/m]
%                 column 2 contains the qc of branch start point [veh/s]
%
% OUTPUTS
%---- q : scalar or vector (same size as k), vehicle flow [veh/s]

scalingfactor = param{1};

kc = param{2}(:,1); % abscissa of branch start point [veh/m]
qc = param{2}(:,2); % ordinate of branch start point [veh/s]

Npts = length(kc);
wc = zeros(Npts-1,1);
qc0 = zeros(Npts-1,1);
for i = 1:(Npts-1)
    wc(i) = (qc(i+1) - qc(i))/(kc(i+1) - kc(i));
    qc0(i) = (qc(i)*kc(i+1) - qc(i+1)*kc(i))/(kc(i+1) - kc(i));
end

wc = wc*ones(1,length(k));
qc0 = qc0*ones(1,length(k));
k = ones(size(wc,1),1)*k;

q = min(wc.*k + scalingfactor.*qc0);

end