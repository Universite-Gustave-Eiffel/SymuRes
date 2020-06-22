function [Texitnew, isupnew] = supplydiscr(Texit,Tlastexit,dN,toutS,qoutS,isup)
% suptime = supplydiscr(Texit,Tlastexit,dN,toutS,qoutS,isup)
% Return the minimum exit time allowed by supply for a given vehicle number
% discretization step dN (dN might either < 1 or > 1)
%
% Nov 2019 - Guilhem Mariotte
%
% INPUTS
%---- Texit     : scalar, desired exit time [s]
%---- Tlastexit : scalar, last exit time [s]
%---- dN        : scalar, vehicle step [veh]
%---- toutS     : vector, times of supply changes [s]
%---- qoutS     : vector, supply values [veh/s]
%---- isup      : integer, current index in the supply values
%
% OUTPUTS
%---- Texitnew : scalar, effective exit time [s]
%---- isupnew  : integer, new current index in the supply values

Nsup = length(qoutS); % number of supply changes

suptime = Tlastexit + dN/qoutS(isup);
Texitnew = max([Texit suptime]);

if isup < Nsup
    dN0 = dN/(Texitnew - Tlastexit)*(toutS(isup+1) - Tlastexit);
    
    while isup < Nsup && toutS(isup+1) < Texitnew
        suptime = toutS(isup+1) + (dN - dN0)/qoutS(isup+1);
        Texitnew = max([Texit suptime]);
        isup = isup + 1;
        if isup < Nsup
            dN0 = dN0 + qoutS(isup)*(toutS(isup+1) - toutS(isup));
        end
    end
end

isupnew = isup;

end