function flow = mergeNd(demand,capacity,control)
% flow = mergeNd(demand,capacity,control)
% Merge function to arbitrate between N incoming flows on a single capacity
% Based on the merge model used in the Meso LWR model by Leclercq & Becarie (2012)
% Extension of the fair merge of Daganzo.
% v00 - exact solution
%
% INPUTS
%---- demand   : vector, demand flows [veh/s]
%---- capacity : scalar, total max allowed flow (total capacity) [veh/s]
%---- control  : vector, each control variable i gives preference to demand i
%                the sum of all variables must be equal to 1
%
% OUTPUTS
%---- flow : vector, same size as demand, return each allowed flow i corresponding to demand i [veh/s]

Nmerge = length(demand); % number of incoming flows

flow = zeros(1,Nmerge);

demand = (0 <= demand).*demand;
capacity = (0 <= capacity).*capacity;
control = (0 <= control).*control;
control = control./sum(control);

eps = 1e-6; % epsilon to test equality in flow [veh/s]

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
            [flowi, imin] = min([demand(i) control(i)/congcontroltot*(capacity - freeflowtot)]);
            flow(i) = flowi;
            if imin == 1
                % The inflow is free-flow
                freeflowtot2 = freeflowtot2 + flow(i);
                Nserved = Nserved + 1;
            else
                % The inflow is congested
                congcontroltot2 = congcontroltot2 + control(i);
                unservedflows2 = [unservedflows2 i];
            end
        else
            flow(i) = 0;
            Nserved = Nserved + 1;
        end
    end
    unservedflows = unservedflows2;
    freeflowtot = freeflowtot + freeflowtot2;
    congcontroltot = congcontroltot2;
    if sum(flow) >= (1 - eps)*capacity
        Nserved = Nmerge;
    end
end

end