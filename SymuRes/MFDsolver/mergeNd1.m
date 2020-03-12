function [inflow, outflow] = mergeNd1(demand,supply)
% [inflow, outflow] = mergeNd1(demand,supply)
% Apply the fair merge of Daganzo in N dimensions with identical merge
% coefficients (equiprobability for each inflow).
%
% INPUTS
%---- demand   : vector, incoming demand flows [veh/s]
%---- supply   : vector, outgoing supply flows [veh/s]
%
% OUTPUTS
%---- inflow  : vector, same size as demand, return each allowed flow i corresponding to demand i [veh/s]
%---- outflow : vector, same size as supply, return each allowed flow i corresponding to supply i [veh/s]

Nmerge = length(demand); % number of incoming flows

control = ones(1,Nmerge)./Nmerge;

inflow = mergeNd(demand,supply,control);
outflow = supply;

end