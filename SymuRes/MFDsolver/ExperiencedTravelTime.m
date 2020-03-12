function TT = ExperiencedTravelTime(t,nin,nout)
% TT = ExperiencedTravelTime(t,nin,nout)
% Calculate the estimated travel time at t for the vehicle exiting at t
% Estimation based on the cumulative entering and exiting vehicle curves
%
% INPUTS
%---- t     : row vector, time [s] 
%---- nin   : row vector, cumulative number of vehicles entering the network, same size as t
%---- nout  : row vector, cumulative number of vehicles exiting the network, same size as t
%
% OUTPUTS
%---- TT : row vector, travel time [s] at t experienced by the vehicle exiting at t

Nt = length(t);
TT = zeros(1,Nt);

for i = 2:Nt
    j = i;
    while j > 2 && nin(j) > nout(i)
        j = j - 1;
    end
    TT(i) = t(i) - (t(j-1) + (t(j) - t(j-1))*(nout(i) - nin(j-1))/(nin(j) - nin(j-1)));
end

end