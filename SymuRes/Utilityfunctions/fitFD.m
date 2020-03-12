function q = fitFD(k,p,Npts,fctFD)
% q = pwlinearFD2(k,p,Npts)
% Modify input arguments of a given FD function
%
% INPUTS
%---- k     : scalar or vector, vehicle density [veh/m], range of k must be between 0 and kj
%---- p     : vector, list of parameters for the FD function
%             p(1) is kj, p(2:Npts+1) are the densities kc, the following values are the flows qc
%---- Npts  : integer, number of middle points to define the FD function
%---- fctFD : Matlab function, FD function of the form fctFD(k,kj,pts)
%
% OUTPUTS
%---- q : scalar or vector (same size as k), vehicle flow [veh/s]

Nparam = 2*Npts + 1;
pts = zeros(Npts,2);
pts(:,1) = p(2:(Npts+1));
pts(:,2) = p((Npts+2):Nparam);

q = fctFD(k,p(1),pts);

end