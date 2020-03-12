function [globalinflow, outflow] = divergeNd(globaldemand,demandratio,supply)
% divergeNd(globaldemand,demandratio,supply)
% Diverge function to arbitrate between 1 incoming flow on M outgoing links
% Model from Tampere et al. (2011)
%
% INPUTS
%---- globaldemand : scalar, global demand into the node [veh/s]
%---- demandratio  : vector, incoming demand turning fractions [-], sum(demandratio) must be 1
%---- supply       : vector, same size as demandratio, outgoing supply flows [veh/s]
%
% OUTPUTS
%---- globalinflow  : scalar, return global inflow in the node [veh/s]
%---- outflow       : vector, same size as supply, return each allowed flow i corresponding to supply i [veh/s]

globalinflow = min([globaldemand supply./demandratio]);
outflow = demandratio.*globalinflow;

end