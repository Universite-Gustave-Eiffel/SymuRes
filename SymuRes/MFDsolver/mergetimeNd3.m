function [time, index] = mergetimeNd3(demandtimes,supplytimes)
% index = mergetimeNd(demandtimes,supplytimes,control)
% Merge function to arbitrate between N entry times with a supply time for
% each entry.
% Select the min effective time = min(demand,supply)
%
% INPUTS
%---- demandtimes : vector, demand entering time per entry [s]
%---- supplytimes : vector, same size, allowed entering time per entry [s]
%
% OUTPUTS
%---- time  : scalar, effective entering time
%---- index : integer, return the entry index for which the supply slot is
%             allowed

indexset = find(supplytimes == -Inf);
supplytimes(indexset) = -ones(1,length(indexset));

% Elmininate entries with infinite supply
indexset = find(supplytimes < Inf);

if isempty(indexset)
    % No entry allowed
    index = 1;
    time = Inf;
else
    % Set of the min effective times (max of demand and supply)
    efftimes = (demandtimes > supplytimes).*demandtimes + (demandtimes <= supplytimes).*supplytimes;
    indexsetmin = find(efftimes(indexset) == min(efftimes(indexset)));
    Nset = length(indexsetmin);
    % Uniform random draw of an entry among the min times
    i = randi(Nset);
    index = indexset(indexsetmin(i));
    time = max([demandtimes(index) supplytimes(index)]);
end

end