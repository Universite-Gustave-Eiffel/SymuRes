function [time, index] = mergetimeNd2(demandtimes,supplytimes)
% index = mergetimeNd(demandtimes,supplytimes,control)
% Merge function to arbitrate between N entry times with a supply time for
% each entry.
% Based on the merge model from Jin & Zhang (2003), pro-rata merge.
%
% INPUTS
%---- demandtimes : vector, demand entering time per entry [s]
%---- supplytimes : vector, same size, allowed entering time per entry [s]
%
% OUTPUTS
%---- time  : scalar, effective entering time
%---- index : integer, return the entry index for which the supply slot is
%             allowed

% Elmininate entries with infinite supply
indexset = find(supplytimes < Inf);

if isempty(indexset)
    % No entry allowed
    index = 1;
    time = Inf;
else
    % Set of the min demand times
    indexsetmin = find(demandtimes(indexset) == min(demandtimes(indexset)));
    Nset = length(indexsetmin);
    % Uniform random draw of an entry among the min times
    i = randi(Nset);
    index = indexset(indexsetmin(i));
    time = max([demandtimes(index) supplytimes(index)]);
end

end