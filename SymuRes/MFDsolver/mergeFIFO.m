function flow = mergeFIFO(itime,t,nincurrent,nindem,demand,capacity)
% flow = mergeFIFO(itime,t,nincurrent,nindem,demand,capacity)
% Merge function to arbitrate between N incoming flows on a single
% capacity. FIFO rule to allocate the capacity to each flow: pro-rata
% splitting rule of the flows when they enter their respective queues.
%
% Nov 2019 - Guilhem Mariotte
%
% INPUTS
%---- itime      : integer, current time index in t
%---- t          : Nt-size vector, simulation time [s]
%---- nincurrent : N-size vector, current effective cumulative counts of
%                  each flow (nb of veh that already transfered) [veh]
%---- nindem     : N-by-Nt matrix, demand cumulative counts of each flow at
%                  each time step [veh]
%---- demand     : N-size vector, demand flows [veh/s]
%---- capacity   : scalar, total max allowed flow (total capacity) [veh/s]
%
% OUTPUTS
%---- flow : N-size vector, return each allowed flow corresponding to each
%            demand [veh/s]

Nflows = size(nindem,1); % number of incoming flows
Nt = length(t);

flow = zeros(1,Nflows);

demand = (0 <= demand).*demand;
capacity = (0 <= capacity).*capacity;
nindemtot = sum(nindem,1);
demandtot = sum(demand);

if itime > 1
    % Total effective inflow (allowed to enter)
    flowtot = min([demandtot capacity]);
    % Total cumulative nb of veh (num of the last entered veh)
    Ntot = sum(nincurrent) + (t(itime) - t(itime-1))*flowtot;
    
    j = itime;
    while j > 2 && nindemtot(j) > Ntot
        j = j - 1;
    end
    % Time at which this veh entered the queue: t0 = t(j) + Dt0
    itime0 = min([j Nt-1]);
    Dt0 = (t(itime0+1) - t(itime0))*(Ntot - nindemtot(itime0))/(nindemtot(itime0+1) - nindemtot(itime0));
    
    for i = 1:Nflows
        % Num of the veh in route i that entered its queue i at t0
        N0 = nindem(i,itime0) + Dt0*(nindem(i,itime0+1) - nindem(i,itime0))/(t(itime0+1) - t(itime0));
        % Allowed inflow for route i based on this num
        flow(i) = (N0 - nincurrent(i))/(t(itime) - t(itime-1));
    end
else
    flow = demand;
end

end