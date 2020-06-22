function departuretime = demdiscr(Tsimul,dN,tinD,qinD)
% departuretime = demdiscr(Tsimul,dN,tinD,qinD)
% Return the departure times for a given vehicle number discretization step
% dN (dN might either < 1 or > 1)
%
% Nov 2019 - Guilhem Mariotte
%
% INPUTS
%---- Tsimul : scalar, simulation duration [s]
%---- dN     : scalar, vehicle step [veh]
%---- tinD   : vector, times of demand changes [s]
%---- qinD   : vector, same size as tinD, demand values [veh/s]
%
% OUTPUTS
%---- departuretime : vector, departure times in chronological order [s]

i = 1;
while i <= length(tinD) && tinD(i) <= Tsimul
    i = i + 1;
end
i = max([1 i-1]);
Ndem = i; % number of demand changes
qinD = qinD(1:i);
tinD = [tinD(1:i) Tsimul];

Ndep = 0; % number of departures
for i = 1:Ndem
    Ndep = Ndep + qinD(i)*(tinD(i+1) - tinD(i))/dN;
end
Ndep = max([1 round(Ndep)]);

departuretime = Inf*ones(1,Ndep);
if qinD(1) > 0
    departuretime(1) = 0;
end

demchange = 0; % token to indicate a change in demand
t0 = 0;
dN0 = 0;
id = 1;
for i = 1:Ndem
    while departuretime(id)+dN/qinD(i) <= tinD(i+1) && demchange == 0
        departuretime(id+1) = departuretime(id) + dN/qinD(i);
        t0 = departuretime(id+1);
        dN0 = 0;
        id = id + 1;
    end
    
    if i < Ndem
        dN0 = dN0 + qinD(i)*(tinD(i+1) - t0);
        departuretime(id+1) = tinD(i+1) + (dN - dN0)/qinD(i+1); % forecasted departure time
        
        if departuretime(id+1) > tinD(i+2) % if a demand change occurs before
            t0 = tinD(i+1);
            demchange = 1;
        else % else t0, dN0 are reset, id is incremented (the previous departure time is valid)
            demchange = 0;
            t0 = departuretime(id+1);
            dN0 = 0;
            id = id + 1;
        end
    end
    
    %     Ni = floor(qinD(i)*(tinD(i+1) - tinD(i))/dN);
    %     departuretime((id+1):(id+Ni)) = departuretime(id) + (1:Ni).*dN/qinD(i);
end

departuretime = departuretime(1:id);

if departuretime(1) == Inf
    if length(departuretime) > 1
        departuretime = departuretime(2:end);
    else
        departuretime = [];
    end
end

end