function [time, index] = selectmintime(demandtimes,supplytimes)
% [time, index] = selectmintime(demandtimes,supplytimes)
% Select the minimum effective time for N demand times with their
% corresponding supply times.
% min eff time = min(max(demandtimes,supplytimes))
% (a uniform random draw is used if the min is not unique)
%
% Nov 2019 - Guilhem Mariotte
%
% INPUTS
%---- demandtimes : N-size vector, demand entering time per entry [s]
%---- supplytimes : N-size vector, allowed entering time per entry [s]
%
% OUTPUTS
%---- time  : scalar, selected effective entry time [s]
%---- index : integer, index of the selected time in the original list

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