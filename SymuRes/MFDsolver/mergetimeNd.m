function [time, index] = mergetimeNd(demandtimes,supplytimes,control)
% index = mergetimeNd(demandtimes,supplytimes,control)
% Merge function to arbitrate between N entry times with a supply time for
% each entry.
% Based on the merge model used in the Meso LWR model by Leclercq & Becarie (2012)
% Extension of the fair merge of Daganzo.
%
% INPUTS
%---- demandtimes : vector, demand entering time per entry [s]
%---- supplytimes : vector, same size, allowed entering time per entry [s]
%---- control     : vector, same size, each control variable i gives
%                   preference to demand i, the sum of all variables must be equal to 1
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
    % Set of congested entries
    congset = find(demandtimes(indexset) < supplytimes(indexset));
    
    if isempty(congset) % all entries in free-flow
        % Set of the min demand times
        indexsetmin = find(demandtimes(indexset) == min(demandtimes(indexset)));
        Nset = length(indexsetmin);
        % Uniform random draw of an entry among the min times
        i = randi(Nset);
        index = indexset(indexsetmin(i));
        time = demandtimes(index);
        
    else % some entries congested
        % Probabilities to select the congested entries
        Nset = length(congset);
        controltot = sum(control(indexset(congset)));
        proba = zeros(1,Nset);
        for i = 1:Nset
            proba(i) = control(indexset(congset(i)))/controltot;
        end
        % Random draw of an entry according to these probabilities
        index = indexset(congset(Nset)); % choose a priori the last index
        draw = rand; % random draw
        minbound = 0;
        for i = 1:(Nset-1)
            if minbound <= draw && draw < minbound+proba(i)
                index = indexset(congset(i)); % may choose another index depending on the random draw
            end
            minbound = minbound + proba(i);
        end
        time = supplytimes(index);
    end
end

end