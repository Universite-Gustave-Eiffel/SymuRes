function [tval, fval] = fstair3(Dt,tval0,fval0)
% [tval, fval] = fstair3(Dt,tval0,fval0)
% Return function values by stairs
%
% INPUTS
%---- Dt    : scalar, width of the stairs
%---- tval0 : vector, time values
%---- fval0 : vector, same size as tval0, function values
%
% OUTPUTS
%---- tval : vector, size Nstairs+1, time values at the stair edges
%---- fval : vector, size Nstairs+1, stair function values

Nt0 = length(tval0);
tstart = tval0(1);
tend = tval0(end);
Nstairs = floor((tend - tstart)/Dt);

tval = zeros(1,Nstairs+1);
fval = zeros(1,Nstairs+1);

tval(1) = tstart;
stairvalue = 0;
i0 = 1;
for i = 1:Nstairs
    tval(i+1) = tval(i) + Dt;
    
    totval = 0;
    Nval = 0;
    while i0 <= Nt0 && tval0(i0) < tval(i+1)
        totval = totval + fval0(i0);
        Nval = Nval + 1;
        i0 = i0 + 1;
    end
    if Nval ~= 0
        stairvalue = totval/Nval;
    end
    fval(i) = stairvalue;
end

fval(Nstairs+1) = stairvalue;

end