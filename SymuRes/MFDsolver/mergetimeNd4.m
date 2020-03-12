function supplytimes = mergetimeNd4(Ltrip,demandtimes,lastentrytimes,capacity,control)
% supplytimes = mergetimeNd4(Ltrip,demandtimes,lastentrytimes,capacity,control)
% Merge function to arbitrate between N entry times with a supply time for
% each entry.
% Based on the merge model used in the Meso LWR model by Leclercq & Becarie (2012)
% Extension of the fair merge of Daganzo.
% Based on the flow formulation of the algorithm, return supply times
%
% INPUTS
%---- Ltrip          : vector, trip lengths corresponding to each entry
%                      (only required for a merge in production, put [1 1...] otherwise)
%---- demandtimes    : vector, demand entry time per entry [s]
%---- lastentrytimes : vector, last entry times per entry [s]
%---- capacity       : scalar, total max allowed flow (total capacity) [veh/s]
%---- control        : vector, same size, each control variable i gives
%                      preference to demand i, the sum of all variables must be equal to 1
%
% OUTPUTS
%---- supplytimes : vector, same size, allowed entry time per entry [s]

Nmerge = length(demandtimes); % number of incoming flows

supplytimes = -Inf*ones(1,Nmerge);
flows = zeros(1,Nmerge);

eps = 1e-3; % epsilon to test equality in flow [veh/s]

unservedflows = 1:Nmerge;
freeflowtot = 0;
congcontroltot = 1;
Nserved = 0;

while Nserved < Nmerge
    unservedflows2 = [];
    freeflowtot2 = 0;
    congcontroltot2 = 0;
    for i = unservedflows
        if control(i) > 0
            supply = control(i)/congcontroltot*(capacity - freeflowtot)/Ltrip(i); % always a flow value
            supplytimes(i) = lastentrytimes(i) + 1/supply;
            if demandtimes(i) > lastentrytimes(i)
                flows(i) = Ltrip(i)*min([1/(demandtimes(i)-lastentrytimes(i)) supply]); % flow or prod value
            else
                flows(i) = Ltrip(i)*supply;
            end
            if demandtimes(i) >= supplytimes(i)
                % The inflow is free-flow
                freeflowtot2 = freeflowtot2 + flows(i);
                Nserved = Nserved + 1;
            else
                % The inflow is congested
                congcontroltot2 = congcontroltot2 + control(i);
                unservedflows2 = [unservedflows2 i];
            end
        else
            supplytimes(i) = Inf;
            Nserved = Nserved + 1;
        end
    end
    unservedflows = unservedflows2;
    freeflowtot = freeflowtot + freeflowtot2; % flow or prod value
    congcontroltot = congcontroltot2;
    if sum(flows) >= capacity - eps
        Nserved = Nmerge;
    end
end

end