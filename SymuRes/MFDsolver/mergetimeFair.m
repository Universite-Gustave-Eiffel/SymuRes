function supplytimes = mergetimeFair(Ltrip,demandtimes,lastentrytimes,dN,capacity,splitcoeff)
% supplytimes = mergetimeFair(Ltrip,demandtimes,lastentrytimes,dN,capacity,splitcoeff)
% Merge function to arbitrate between N entry demand times on a single
% capacity. Extension of the fair merge of Daganzo (1995), merge model used
% in the Meso LWR model by Leclercq & Becarie (2012). Based on the flow
% formulation of the algorithm, return supply times.
%
% Nov 2019 - Guilhem Mariotte
%
% INPUTS
%---- Ltrip          : N-size vector, trip lengths corresponding to each entry
%                      (only required for a merge in production, put [1 1...] otherwise)
%---- demandtimes    : N-size vector, demand entry time per entry [s]
%---- lastentrytimes : N-size vector, last entry times per entry [s]
%---- dN             : integer > 0, vehicle increment, generally of 1 [veh]
%---- capacity       : scalar, total max allowed flow (total capacity) [veh/s]
%---- splitcoeff     : N-size vector, coefficients that define the splitting
%                      rule among the demand flows (sum coeffs = 1)
%
% OUTPUTS
%---- supplytimes : N-size vector, allowed entry time per entry [s]

Nflows = length(demandtimes); % number of incoming flows

supplytimes = -Inf*ones(1,Nflows);
flows = zeros(1,Nflows);

eps = 1e-3; % epsilon to test equality in flow [veh/s]

unservedflows = 1:Nflows;
freeflowtot = 0;
congcoefftot = 1;
Nserved = 0;

while Nserved < Nflows
    unservedflows2 = [];
    freeflowtot2 = 0;
    congcoefftot2 = 0;
    for i = unservedflows
        if splitcoeff(i) > 0
            supply = splitcoeff(i)/congcoefftot*(capacity - freeflowtot)/Ltrip(i); % always a flow value
            supplytimes(i) = lastentrytimes(i) + dN/supply;
            if demandtimes(i) > lastentrytimes(i)
                flows(i) = Ltrip(i)*min([dN/(demandtimes(i)-lastentrytimes(i)) supply]); % flow or prod value
            else
                flows(i) = Ltrip(i)*supply;
            end
            if demandtimes(i) >= supplytimes(i)
                % The inflow is free-flow
                freeflowtot2 = freeflowtot2 + flows(i);
                Nserved = Nserved + 1;
            else
                % The inflow is congested
                congcoefftot2 = congcoefftot2 + splitcoeff(i);
                unservedflows2 = [unservedflows2 i];
            end
        else
            supplytimes(i) = Inf;
            Nserved = Nserved + 1;
        end
    end
    unservedflows = unservedflows2;
    freeflowtot = freeflowtot + freeflowtot2; % flow or prod value
    congcoefftot = congcoefftot2;
    if sum(flows) >= capacity - eps
        Nserved = Nflows;
    end
end

supplytimes = max([demandtimes; supplytimes]);

end