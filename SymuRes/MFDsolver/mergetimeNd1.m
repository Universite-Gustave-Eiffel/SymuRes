function [time, index] = mergetimeNd1(demandtimes,supplytimes)
% [time, index] = mergetimeNd1(demandtimes,supplytimes)
% Apply the fair merge of Daganzo in N dimensions with identical merge
% coefficients (equiprobability for each inflow).
%
% INPUTS
%---- demandtimes : vector, demand entering time per entry [s]
%---- supplytimes : vector, same size, allowed entering time per entry [s]
%
% OUTPUTS
%---- time  : scalar, effective entering time
%---- index : integer, return the entry index for which the supply slot is
%             allowed

Nmerge = length(demandtimes); % number of incoming flows

control = ones(1,Nmerge)./Nmerge;

[time, index] = mergetimeNd(demandtimes,supplytimes,control);

end