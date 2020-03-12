function P = MFDest(n,param,MFDfct)
% P = MFDest(n,param,MFDfct)
% Add constraints on the MFD function model for MFD fit. Mainly to be used
% with parabolic function models.
%
% INPUTS
%---- n      : vector, values of accumulation [veh]
%---- param  : vector, parameters for the MFD function model, must begin
%              with [nj nc ...]
%---- MFDfct : function handle, prod-MFD function model P(n,parameters)
%
% OUTPUTS
%---- P : vector, same size as n, values of production [veh.m/s]

nj = param(1);
nc = param(2);
if ~(2*nc <= nj && nj < 3*nc)
    P = zeros(1,length(n));
else
    P = MFDfct(n,param);
end

end