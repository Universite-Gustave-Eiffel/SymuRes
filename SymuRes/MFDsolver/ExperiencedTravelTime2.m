function TT = ExperiencedTravelTime2(t,nin,nout,it1,it2)
% TT = ExperiencedTravelTime2(t,nin,nout,it1,it2)
% Calculate the estimated travel time at t for the vehicle exiting at t
% over the period [t(it1) t(it2)].
% Estimation based on the cumulative entering and exiting vehicle curves
%
% INPUTS
%---- t     : row vector, time [s] 
%---- nin   : row vector, cumulative number of vehicles entering the network, same size as t
%---- nout  : row vector, cumulative number of vehicles exiting the network, same size as t
%---- it1   : integer, indice of time t, start of the period
%---- it2   : integer, indice of time t, end of the period
%
% OUTPUTS
%---- TT : row vector, travel time [s] at t experienced by the vehicle exiting at t

Nt = length(t);

% default period = total time frame, if not mentioned
if nargin == 3
    it1 = 1;
    it2 = Nt;
end

NTT = it2 - it1 + 1;
TT = zeros(1,NTT);

for i = 2:NTT
    j = i + it1 - 1;
    while j > 2 && nin(j) > nout(i+it1-1)
        j = j - 1;
    end
    TT(i) = t(i+it1-1) - (t(j-1) + (t(j) - t(j-1))*(nout(i+it1-1) - nin(j-1))/(nin(j) - nin(j-1)));
end

if it1 > 1
    i = 1;
    j = i + it1 - 1;
    while j > 2 && nin(j) > nout(i+it1-1)
        j = j - 1;
    end
    TT(i) = t(i+it1-1) - (t(j-1) + (t(j) - t(j-1))*(nout(i+it1-1) - nin(j-1))/(nin(j) - nin(j-1)));
end

end