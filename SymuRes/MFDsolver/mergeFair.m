function flow = mergeFair(demand,capacity,splitcoeff)
% flow = mergeFair(demand,capacity,splitcoeff)
% Merge function to arbitrate between N incoming flows on a single
% capacity. Extension of the fair merge of Daganzo (1995), merge model used
% in the Meso LWR model by Leclercq & Becarie (2012).
%
% Nov 2019 - Guilhem Mariotte
%
% INPUTS
%---- demand     : N-size vector, demand flows [veh/s]
%---- capacity   : scalar, total max allowed flow (total capacity) [veh/s]
%---- splitcoeff : N-size vector, coefficients that define the splitting
%                  rule among the demand flows (sum coeffs = 1)
%
% OUTPUTS
%---- flow : N-size vector, return each allowed flow corresponding to each
%            demand [veh/s]

Nflows = length(demand); % number of incoming flows

flow = zeros(1,Nflows);

demand = (0 <= demand).*demand;
capacity = (0 <= capacity).*capacity;
splitcoeff = (0 <= splitcoeff).*splitcoeff;
splitcoeff = splitcoeff./sum(splitcoeff);

eps = 1e-6; % epsilon to test equality in flow [veh/s]

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
            [flowi, imin] = min([demand(i) splitcoeff(i)/congcoefftot*(capacity - freeflowtot)]);
            flow(i) = flowi;
            if imin == 1
                % The inflow is free-flow
                freeflowtot2 = freeflowtot2 + flow(i);
                Nserved = Nserved + 1;
            else
                % The inflow is congested
                congcoefftot2 = congcoefftot2 + splitcoeff(i);
                unservedflows2 = [unservedflows2 i];
            end
        else
            flow(i) = 0;
            Nserved = Nserved + 1;
        end
    end
    unservedflows = unservedflows2;
    freeflowtot = freeflowtot + freeflowtot2;
    congcoefftot = congcoefftot2;
    if sum(flow) >= (1 - eps)*capacity
        Nserved = Nflows;
    end
end

end