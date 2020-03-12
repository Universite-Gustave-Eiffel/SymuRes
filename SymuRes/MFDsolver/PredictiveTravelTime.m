function TT = PredictiveTravelTime(t,nin,nout)
% TT = PredictiveTravelTime(t,nin,nout)
% Calculate the travel time at t for the vehicle entering at t
% Calculation based on the cumulative entering and exiting vehicle curves
%
% INPUTS
%---- t    : row vector, time [s] 
%---- nin  : row vector, cumulative number of vehicles entering the network, same size as t
%---- nout : row vector, cumulative number of vehicles exiting the network, same size as t
%
% OUTPUTS
%---- TT : row vector, travel time [s] at t experienced by the vehicle entering at t

Nt = length(t);
TT = zeros(1,Nt);

for i = 1:(Nt-1) % loop on all times
    j = i + 1;
    while j < Nt && nin(i) > nout(j)
        j = j + 1;
    end
    TT(i) = (t(j-1) + (t(j) - t(j-1))*(nin(i) - nout(j-1))/(nout(j) - nout(j-1))) - t(i);
end

end