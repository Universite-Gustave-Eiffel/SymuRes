function [inflow, outflow] = mergeNd2(demand,supply)
% [inflow, outflow] = mergeNd2(demand,supply)
% Merge function to arbitrate between N incoming flows on M outgoing links
% Model from Jin & Zhang (2003), pro-rata merge.
% No routes are considered for outflows, inflows and outflows are
% determined proportionally to demand and supply.
%
% INPUTS
%---- demand   : vector, incoming demand flows [veh/s]
%---- supply   : vector, outgoing supply flows [veh/s]
%
% OUTPUTS
%---- inflow  : vector, same size as demand, return each allowed flow i corresponding to demand i [veh/s]
%---- outflow : vector, same size as supply, return each allowed flow i corresponding to supply i [veh/s]

globaldemand = sum(demand);
globalsupply = sum(supply);
globalflow = min([globaldemand globalsupply]);

if globaldemand > 0
    inflow = demand./globaldemand.*globalflow;
else
    inflow = 0.*demand;
end

if globalsupply > 0
    outflow = supply./globalsupply.*globalflow;
else
    outflow = 0.*supply;
end

end