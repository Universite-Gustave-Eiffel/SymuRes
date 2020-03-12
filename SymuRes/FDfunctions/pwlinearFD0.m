function q = pwlinearFD0(k,p,Npts)
% q = pwlinearFD2(k,p,Npts)
% Compute the piecewise linear Fundamental Diagram function 
%
% INPUTS
%---- k    : scalar or vector, vehicle density [veh/m], range of k must be between 0 and kj
%---- p    : vector, list of parameters for the piecewise linear FD
%            p(1) is kj, p(2:Npts+1) are the densities kc, the following values are the flows qc
%---- Npts : integer, number of points to define the piecewise linear FD
%
% OUTPUTS
%---- q : scalar or vector (same size as k), vehicle flow [veh/s]

Nparam = 2*Npts + 1;
pts = zeros(Npts,2);
pts(:,1) = p(2:(Npts+1));
pts(:,2) = p((Npts+2):Nparam);

q = pwlinearFD(k,[pts; p(1) 0]);

end