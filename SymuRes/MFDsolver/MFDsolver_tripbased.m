%% MFD SOLVER
%--------------------------------------------------------------------------
% Multi-reservoir MFD-based traffic flow solver
% Trip-based model
%
% Nov 2019 - Guilhem Mariotte
%
% References:
% Mariotte et al. (TR part B, 2017)
% Mariotte & Leclercq (TR part B, 2019)
% Mariotte et al. (TR part B, 2020)


%% Vehicle creation
%--------------------------------------------------------------------------

% Simulation attributes
SimulationDuration = Simulation.Duration;
TimeStep = Simulation.TimeStep;
MFDfct = Simulation.MFDfct;
Entryfct = Simulation.Entryfct;

% Assignment period time window
Temp_StartTimeID = floor(Assignment.CurrentTime/TimeStep) + 1;
Temp_EndTimeID = min([floor(Assignment.Periods(Assignment.CurrentPeriodID+1)/TimeStep) NumTimes-1]);
Temp_numtimesperiod = Temp_EndTimeID - Temp_StartTimeID + 1;
Temp_StartTime = Assignment.CurrentTime;
Temp_EndTime = min([Assignment.Periods(Assignment.CurrentPeriodID+1) Simulation.Duration]);

% Load from previous assignment period
if Assignment.CurrentPeriodID > 1
    Reservoir = Snapshot.Reservoir;
    Global = Snapshot.Global;
    Vehicle = Snapshot.Vehicle;
    NextEvent = Snapshot.NextEvent;
end

% Append entry times to the routes
Temp_routetimeID = ones(1,NumRoutes);

for iroute = 1:NumRoutes
    od = Route(iroute).ODmacroID;
    
    % Initialize for the first assignment period
    if Assignment.CurrentPeriodID == 1
        Route(iroute).EntryTimes = [];
        Route(iroute).EntryPurposes = {};
    end
    
    if Route(iroute).AssignCoeff > 0
        Temp_istart = length(Route(iroute).EntryTimes) + 1;
        Temp_routetimeID(iroute) = Temp_istart;
        
        for j = 1:length(ODmacro(od).Demand)
            Temp_demtimes = Simulation.Time;
            Temp_demdata = Simulation.TripbasedSimuFactor.*Route(iroute).Demand;
            it1 = 1;
            while it1 <= length(Temp_demtimes) && Temp_demtimes(it1) < Temp_StartTime
                it1 = it1 + 1;
            end
            it1 = max([1 it1-1]);
            it2 = it1;
            while it2 <= length(Temp_demtimes) && Temp_demtimes(it2) < Temp_EndTime
                it2 = it2 + 1;
            end
            it2 = max([it1 it2-1]);
            
            if Assignment.CurrentPeriodID == 1
                Temp_demtimes12 = Temp_demtimes(it1:it2);
                Temp_demdata12 = Temp_demdata(it1:it2);
                Temp_times = demdiscr(Temp_EndTime,1,Temp_demtimes12,Temp_demdata12);
                if isempty(Temp_times)
                    Route(iroute).PrevDemandTime = Temp_demtimes12;
                    Route(iroute).PrevDemandData = Temp_demdata12;
                else
                    Route(iroute).PrevDemandTime = Temp_times(end);
                    Route(iroute).PrevDemandData = Temp_demdata12(end);
                end
            else
                if it1 == it2
                    Temp_demtimes12 = [Route(iroute).PrevDemandTime Temp_StartTime];
                    Temp_demdata12 = [Route(iroute).PrevDemandData Temp_demdata(it1)];
                else
                    Temp_demtimes12 = [Route(iroute).PrevDemandTime Temp_StartTime Temp_demtimes((it1+1):it2)];
                    Temp_demdata12 = [Route(iroute).PrevDemandData Temp_demdata(it1) Temp_demdata((it1+1):it2)];
                end
                Temp_times = demdiscr(Temp_EndTime-Temp_demtimes12(1),1,[0 Temp_demtimes12-Temp_demtimes12(1)],[0 Temp_demdata12]);
                Temp_times = Temp_times + Temp_demtimes12(1);
                if isempty(Temp_times)
                    Route(iroute).PrevDemandTime = Temp_demtimes12;
                    Route(iroute).PrevDemandData = Temp_demdata12;
                else
                    Route(iroute).PrevDemandTime = Temp_times(end);
                    Route(iroute).PrevDemandData = Temp_demdata12(end);
                end
            end
            
            Route(iroute).EntryTimes = [Route(iroute).EntryTimes Temp_times];
            for k = 1:length(Temp_times)
                Route(iroute).EntryPurposes = [Route(iroute).EntryPurposes {ODmacro(od).Demand(j).Purpose}];
            end
            
        end
        
        Temp_iend = length(Route(iroute).EntryTimes);
        Temp_Ntimes = Temp_iend - Temp_istart + 1;
        
        % Sort by ascending entry times
        [Temp_times, Temp_sortindex] = sort(Route(iroute).EntryTimes(Temp_istart:Temp_iend));
        Route(iroute).EntryTimes(Temp_istart:Temp_iend) = Temp_times;
        Route(iroute).NumEntryTimes = length(Route(iroute).EntryTimes);
        
        Temp_purposes = cell(1,Temp_Ntimes);
        for j = 1:Temp_Ntimes
            Temp_purposes{j} = Route(iroute).EntryPurposes{Temp_istart+Temp_sortindex(j)-1};
        end
        Route(iroute).EntryPurposes(Temp_istart:Temp_iend) = Temp_purposes;
    end
end

% keyboard

% Initialize for the first assignment period
if Assignment.CurrentPeriodID == 1
    Global.EntryTimes = [];
    Global.EntryRoutes = [];
    Global.EntryPurposes = {};
    Global.SimulTime = [];
    Global.VehID = [];
end

% Append entry times to the Global structure
Temp_istart = length(Global.EntryTimes) + 1;

for iroute = 1:NumRoutes
    if Route(iroute).AssignCoeff > 0
        Temp_times = Route(iroute).EntryTimes(Temp_routetimeID(iroute):end);
        Temp_purposes = Route(iroute).EntryPurposes(Temp_routetimeID(iroute):end);
        
        Global.EntryTimes = [Global.EntryTimes Temp_times];
        Global.EntryRoutes = [Global.EntryRoutes iroute*ones(1,length(Temp_times))];
        Global.EntryPurposes = [Global.EntryPurposes Temp_purposes];
    end
end

Temp_iend = length(Global.EntryTimes);
Temp_Ntimes = Temp_iend - Temp_istart + 1;

% Sort by ascending entry times
[Temp_times, Temp_sortindex] = sort(Global.EntryTimes(Temp_istart:Temp_iend));
Global.EntryTimes(Temp_istart:Temp_iend) = Temp_times;
Global.NumEntryTimes = length(Global.EntryTimes);

Temp_routes = zeros(1,Temp_Ntimes);
Temp_purposes = cell(1,Temp_Ntimes);
for j = 1:Temp_Ntimes
    Temp_routes(j) = Global.EntryRoutes(Temp_istart+Temp_sortindex(j)-1);
    Temp_purposes{j} = Global.EntryPurposes{Temp_istart+Temp_sortindex(j)-1};
end
Global.EntryRoutes(Temp_istart:Temp_iend) = Temp_routes;
Global.EntryPurposes(Temp_istart:Temp_iend) = Temp_purposes;




% Initialize for the first assignment period
if Assignment.CurrentPeriodID == 1
    NumVeh = Global.NumEntryTimes;
    Vehicle = struct('RouteID',cell(1,NumVeh));
    iveh = 1;
    
    for r = 1:NumRes
        Reservoir(r).DemandEntryTimePerRoute = cell(1,length(Reservoir(r).RoutesID));
        Reservoir(r).DemandEntryVehPerRoute  = cell(1,length(Reservoir(r).RoutesID));
        Reservoir(r).DemandEntryTime = [];
        Reservoir(r).DemandEntryVeh  = [];
    end
else
    
    iveh = length(Vehicle) + 1;
end

% Vehicle creation and assignment on the routes
for i = Temp_istart:Temp_iend % loop on entry events of the current assignment period
    iroute = Global.EntryRoutes(i);
    r = Route(iroute).ResOriginID;
    i_r = Route(iroute).ResRouteIndex(r);
    Temp_Npath = length(Route(iroute).ResPath);
    
    % Vehicle creation and assignment of a trip length
    % corresponding to the chosen route
    Vehicle(iveh).RouteID = iroute;
    Vehicle(iveh).CreationTime = Global.EntryTimes(i);
    Vehicle(iveh).EntryTimes = Inf*ones(1,Temp_Npath);
    Vehicle(iveh).ExitTimes = Inf*ones(1,Temp_Npath);
    Vehicle(iveh).TripLength = zeros(1,Temp_Npath);
    Vehicle(iveh).TraveledDistance = zeros(1,Temp_Npath);
    Vehicle(iveh).WaitingTimes = zeros(1,Temp_Npath);
    Vehicle(iveh).Purpose = Global.EntryPurposes{i};
    Vehicle(iveh).CurrentResID = r;
    Vehicle(iveh).PathIndex = 1;
    Reservoir(r).DemandEntryTimePerRoute{i_r} = [Reservoir(r).DemandEntryTimePerRoute{i_r} Global.EntryTimes(i)];
    Reservoir(r).DemandEntryVehPerRoute{i_r} = [Reservoir(r).DemandEntryVehPerRoute{i_r} iveh];
    Reservoir(r).DemandEntryTime = [Reservoir(r).DemandEntryTime Global.EntryTimes(i)];
    Reservoir(r).DemandEntryVeh = [Reservoir(r).DemandEntryVeh iveh];
    
    iveh = iveh + 1;
end

NumVeh = iveh - 1;


%% Initialization
%--------------------------------------------------------------------------

% A priori max number of events (one entry and one exit by vehicle by
% reservoir)
Global.NumMaxEvents = 2*NumRes*NumVeh;

% Simulation time
Global.SimulTime = [Global.SimulTime zeros(1,2*NumRes*Temp_Ntimes)];
Global.VehID = [Global.VehID zeros(1,2*NumRes*Temp_Ntimes)];

% Initialize for the first assignment period
if Assignment.CurrentPeriodID == 1
    
    for r = 1:NumRes
        Temp_Nroutes = length(Reservoir(r).RoutesID);
        
        % Simulation variable initialization
        if Temp_Nroutes > 0
            Reservoir(r).VehList = [];
            Reservoir(r).FirstVehPerRoute = zeros(1,Temp_Nroutes);
            Reservoir(r).DemandTimeIndexPerRoute = ones(1,Temp_Nroutes);
            Reservoir(r).DemandTimeIndex = 1;
            
            Reservoir(r).LastEntryTime = -Inf;
            Reservoir(r).LastExitTime = -Inf;
            Reservoir(r).LastEntryTimePerRoute = -Inf*ones(1,Temp_Nroutes);
            Reservoir(r).PrevLastEntryTimePerRoute = -Inf*ones(1,Temp_Nroutes);
            Reservoir(r).LastExitTimePerRoute = -Inf*ones(1,Temp_Nroutes);
            Reservoir(r).NextEntryTime = Inf;
            Reservoir(r).NextExitTime = Inf;
            Reservoir(r).NextEntryVehID = 0;
            Reservoir(r).NextExitVehID = 0;
            
            Reservoir(r).SupplyIndex = 1; % index for the current value of supply
            
            Reservoir(r).DesiredEntryTimePerRoute = Inf*ones(1,Temp_Nroutes);
            Reservoir(r).DesiredExitTimePerRoute = Inf*ones(1,Temp_Nroutes);
            Reservoir(r).DesiredEntryTime = Inf;
            Reservoir(r).DesiredExitTime = Inf;
            Reservoir(r).DesiredEntryVehPerRoute = zeros(1,Temp_Nroutes);
            Reservoir(r).DesiredExitVehPerRoute = zeros(1,Temp_Nroutes);
            Reservoir(r).DesiredEntryVeh = 0;
            Reservoir(r).DesiredExitVeh = 0;
            
            Reservoir(r).EntrySupplyTimePerRoute = zeros(1,Temp_Nroutes);
            Reservoir(r).ExitSupplyTimePerRoute = zeros(1,Temp_Nroutes);
            Reservoir(r).EntrySupplyTime = 0;
            Reservoir(r).ExitSupplyTime = 0;
            
            Reservoir(r).EntryTimes = [];
            Reservoir(r).ExitTimes = [];
            Reservoir(r).EntryTimesPerRoute = cell(1,Temp_Nroutes);
            Reservoir(r).ExitTimesPerRoute = cell(1,Temp_Nroutes);
            
            Reservoir(r).CurrentAcc = 0;
            Reservoir(r).CurrentAccPerRoute = zeros(1,Temp_Nroutes);
            Reservoir(r).CurrentNinPerRoute = zeros(1,Temp_Nroutes);
            Reservoir(r).CurrentNoutPerRoute = zeros(1,Temp_Nroutes);
            Reservoir(r).CurrentMeanSpeed = Reservoir(r).FreeflowSpeed;
            Reservoir(r).MeanSpeed = zeros(1,NumTimes);
            Reservoir(r).InternalProd = 0;
            Reservoir(r).EntryNodeSupplyPerRoute = Inf*ones(1,Temp_Nroutes);
            
            Reservoir(r).MergeCoeffPerRoute = ones(1,Temp_Nroutes);
            Reservoir(r).ExitCoeffPerRoute = ones(1,Temp_Nroutes);
        end
    end
    
    % Origin and destination reservoir lists
    Temp_OriResList = [];
    Temp_DestResList = [];
    for iroute = 1:NumRoutes
        Route(iroute).TravelTime = zeros(1,NumTimes);
        
        if Route(iroute).AssignCoeff > 0 % if the route is used
            o = Route(iroute).ResOriginID;
            d = Route(iroute).ResDestinationID;
            Temp_OriResList = [Temp_OriResList o];
            Temp_DestResList = [Temp_DestResList d];
        end
    end
    Global.OriResList = unique(Temp_OriResList);
    Global.DestResList = unique(Temp_DestResList);
    
    for inode = 1:NumMacroNodes
        MacroNode(inode).Capacity.TimeIndex = 1;
        MacroNode(inode).Supply = MacroNode(inode).Capacity.Data(1);
    end
    
    iveh = 1;
    Global.CurrentTimeID = 1;
    NextEvent.Time = Vehicle(iveh).CreationTime;
    NextEvent.ElapsedTime = 0;
    NextEvent.Type = 1; % 1: entry (vehicle creation) / 2: exit
    NextEvent.VehID = iveh; % current veh ID for the considered event
    NextEvent.ResList = []; % reservoir IDs where accumulation and thus mean speed just changed
    
end

% Initialize at each period
for r = 1:NumRes
    Reservoir(r).MeanSpeed2 = zeros(1,Global.NumMaxEvents);
    Reservoir(r).InflowSupply = zeros(1,Global.NumMaxEvents);
    Reservoir(r).EntryTripLength = zeros(1,Global.NumMaxEvents);
    Reservoir(r).ProdSupply = zeros(1,Global.NumMaxEvents);
end
for iroute = 1:NumRoutes
    Route(iroute).TravelTime2 = zeros(1,Global.NumMaxEvents);
end

% Initialize the desired entry time and next entry in origin reservoirs
for r = Global.OriResList
    % Internal demand
    for i_r = [Reservoir(r).OriRoutesIndex Reservoir(r).EntryRoutesIndex]
        Temp_index = Reservoir(r).DemandTimeIndexPerRoute(i_r);
        if Temp_index <= length(Reservoir(r).DemandEntryTimePerRoute{i_r})
            Reservoir(r).DesiredEntryTimePerRoute(i_r) = Reservoir(r).DemandEntryTimePerRoute{i_r}(Temp_index);
            Reservoir(r).DesiredEntryVehPerRoute(i_r) = Reservoir(r).DemandEntryVehPerRoute{i_r}(Temp_index);
        else
            Reservoir(r).DesiredEntryTimePerRoute(i_r) = Inf;
            Reservoir(r).DesiredEntryVehPerRoute(i_r) = 1;
        end
    end
    Temp_index = Reservoir(r).DemandTimeIndex;
    if Temp_index <= length(Reservoir(r).DemandEntryTime)
        Reservoir(r).DesiredEntryTime = Reservoir(r).DemandEntryTime(Temp_index);
        Reservoir(r).DesiredEntryVeh = Reservoir(r).DemandEntryVeh(Temp_index);
    else
        Reservoir(r).DesiredEntryTime = Inf;
        Reservoir(r).DesiredEntryVeh = 1;
    end
    
    % Next entry in the origin reservoir (beginning of a route)
    Temp_indexes = [];
    for i_r = 1:length(Reservoir(r).RoutesID)
        iroute = Reservoir(r).RoutesID(i_r);
        if r == Route(iroute).ResOriginID
            Temp_indexes = [Temp_indexes i_r];
        end
    end
    if ~isempty(Temp_indexes)
        if strcmp(Simulation.MergeModel,'demfifo') && length(Reservoir(r).EntryNodesIndex) == 1
            % FIFO discipline merge only if one external entry node
            Temp_demandtime = Reservoir(r).DesiredEntryTime; % not necessarily route beginning!
            Temp_supplytime = Reservoir(r).EntrySupplyTime;
            iveh = Reservoir(r).DesiredEntryVeh;
            iroute = Vehicle(iveh).RouteID;
            i_r = Route(iroute).ResRouteIndex(r);
            Temp_time2 = Reservoir(r).LastEntryTime + 1/Reservoir(r).EntryNodeSupplyPerRoute(i_r);
            Reservoir(r).NextEntryTime = max([Temp_demandtime Temp_supplytime Temp_time2]);
            Reservoir(r).NextEntryVehID = iveh;
        else
            % Other merge models
            Temp_demandtimes = Reservoir(r).DesiredEntryTimePerRoute(Temp_indexes);
            Temp_supplytimes = Reservoir(r).EntrySupplyTimePerRoute(Temp_indexes);
            [Temp_time, i_r] = selectmintime(Temp_demandtimes,Temp_supplytimes);
            iveh = Reservoir(r).DesiredEntryVehPerRoute(Temp_indexes(i_r));
            Temp_time2 = Reservoir(r).LastEntryTime + 1/Reservoir(r).EntryNodeSupplyPerRoute(Temp_indexes(i_r));
            Reservoir(r).NextEntryTime = max([Temp_time Temp_time2 CurrentTime]);
            Reservoir(r).NextEntryVehID = iveh;
        end
    else
        Reservoir(r).NextEntryTime = Inf;
        Reservoir(r).NextEntryVehID = 1;
    end
end

eps0 = 0.1;
eps1 = 0.01; % for speed [m/s]
itime = Global.CurrentTimeID;
CurrentTime = NextEvent.Time;
ElapsedTime = NextEvent.ElapsedTime;



%% Simulation loop
%--------------------------------------------------------------------------

while CurrentTime < Temp_EndTime
    
    % Print simulation state
    if floor(NextEvent.VehID/1000) == NextEvent.VehID/1000
        fprintf('%s%3.3f \t %s%i \t %s%i \n','time=',CurrentTime,'nextevent=',NextEvent.Type,'vehID=',NextEvent.VehID)
    end
    
    
    % Update traveled distances
    %--------------------------
    for r = 1:NumRes
        Temp_dist = ElapsedTime*Reservoir(r).CurrentMeanSpeed;
        Temp_Nveh = length(Reservoir(r).VehList);
        for i = 1:Temp_Nveh
            iveh = Reservoir(r).VehList(i);
            i_p = Vehicle(iveh).PathIndex;
            Vehicle(iveh).TraveledDistance(i_p) = Vehicle(iveh).TraveledDistance(i_p) + Temp_dist;
            if Vehicle(iveh).TraveledDistance(i_p) >= Vehicle(iveh).TripLength(i_p)
                Vehicle(iveh).WaitingTimes(i_p) = Vehicle(iveh).WaitingTimes(i_p) + ElapsedTime;
            end
        end
    end
    
    % Simulation time
    Global.SimulTime(itime) = CurrentTime;
    Global.VehID(itime) = NextEvent.VehID;
    
    % Current vehicle, route and reservoir
    iveh = NextEvent.VehID;
    ires = Vehicle(iveh).CurrentResID;
    iroute = Vehicle(iveh).RouteID;
    i_r = Route(iroute).ResRouteIndex(ires);
    Temp_pathindex = Vehicle(iveh).PathIndex;
    NextEvent.ResList = [];
    
    % Current node supply
    for inode = 1:NumMacroNodes
        Temp_index = MacroNode(inode).Capacity.TimeIndex;
        if Temp_index < length(MacroNode(inode).Capacity.Time) && CurrentTime >= MacroNode(inode).Capacity.Time(Temp_index+1)
            Temp_index = Temp_index + 1;
            MacroNode(inode).Capacity.TimeIndex = Temp_index;
            MacroNode(inode).Supply = MacroNode(inode).Capacity.Data(Temp_index);
        end
    end
    
    
    if Vehicle(iveh).EntryTimes(1) < Inf % if not at the beginning of the route
        % Exit of the vehicle
        %--------------------
        
        % Vehicle exit time
        Vehicle(iveh).ExitTimes(Temp_pathindex) = CurrentTime;
        Reservoir(ires).LastExitTime = CurrentTime;
        Reservoir(ires).LastExitTimePerRoute(i_r) = CurrentTime;
        
        % Update reservoir accumulation
        Reservoir(ires).CurrentAcc = Reservoir(ires).CurrentAcc - 1;
        Reservoir(ires).CurrentAccPerRoute(i_r) = Reservoir(ires).CurrentAccPerRoute(i_r) - 1;
        NextEvent.ResList = [NextEvent.ResList ires];
        
        % Update exiting curve
        Reservoir(ires).ExitTimes = [Reservoir(ires).ExitTimes CurrentTime];
        Reservoir(ires).ExitTimesPerRoute{i_r} = [Reservoir(ires).ExitTimesPerRoute{i_r} CurrentTime];
        Reservoir(ires).CurrentNoutPerRoute(i_r) = Reservoir(ires).CurrentNoutPerRoute(i_r) + 1;
        
        % Update reservoir mean speed
        Temp_nr = Reservoir(ires).CurrentAcc;
        Temp_param = Reservoir(ires).MFDfctParam;
        Temp_scf = Simulation.TripbasedSimuFactor;
        if Temp_nr == 0
            Temp_Vr = Reservoir(ires).FreeflowSpeed;
        else
            Temp_Vr = MFDfct(Temp_nr,Temp_scf*Temp_param)/Temp_nr;
        end
        Reservoir(ires).CurrentMeanSpeed = Temp_Vr;
        Reservoir(ires).MeanSpeed2(itime) = Temp_Vr;
        
        % Remove from the global waiting list
        i = 1;
        while i <= length(Reservoir(ires).VehList) && Reservoir(ires).VehList(i) ~= iveh
            i = i + 1;
        end
        Reservoir(ires).VehList(i) = [];
        
        % Find the next first vehicle to exit for the current route
        i = 1;
        while i <= length(Reservoir(ires).VehList) && Vehicle(Reservoir(ires).VehList(i)).RouteID ~= iroute
            i = i + 1;
        end
        if i <= length(Reservoir(ires).VehList)
            Reservoir(ires).FirstVehPerRoute(i_r) = Reservoir(ires).VehList(i);
        else
            Reservoir(ires).FirstVehPerRoute(i_r) = 0; % no next veh found on this route
        end
        
    end
    
    
    if Vehicle(iveh).EntryTimes(1) == Inf || ires ~= Route(iroute).ResDestinationID % if at the beginning of the route or not at destination
        % Entry of the vehicle
        %---------------------
        
        if Vehicle(iveh).EntryTimes(1) < Inf && ires ~= Route(iroute).ResDestinationID % if not at the beginning of the route and not at destination
            % Entry in the next reservoir of the route
            Temp_pathindex = Temp_pathindex + 1;
            ires0 = ires; % previous reservoir ID
            ires = Route(iroute).ResPath(Temp_pathindex); % current reservoir ID
            i_r = Route(iroute).ResRouteIndex(ires);
            Vehicle(iveh).CurrentResID = ires;
            Vehicle(iveh).PathIndex = Temp_pathindex;
        end
        
        % Vehicle entry time
        Vehicle(iveh).EntryTimes(Temp_pathindex) = CurrentTime;
        Reservoir(ires).LastEntryTime = CurrentTime;
        Reservoir(ires).PrevLastEntryTimePerRoute(i_r) = Reservoir(ires).LastEntryTimePerRoute(i_r);
        Reservoir(ires).LastEntryTimePerRoute(i_r) = CurrentTime;
        
        % Set trip length for the current reservoir
        Temp_Ltrip = Reservoir(ires).TripLengthPerRoute(i_r);
        %Temp_LtripStd = Reservoir(ires).TripLengthStdPerODPerPath(o,d,p);
        %Temp_LtripStd = 0.05*Temp_Ltrip;
        %Vehicle(iveh).TripLength(Temp_pathindex) = Temp_Ltrip + Temp_LtripStd*randn;
        Vehicle(iveh).TripLength(Temp_pathindex) = Temp_Ltrip;
        
        % Update reservoir accumulation
        Reservoir(ires).CurrentAcc = Reservoir(ires).CurrentAcc + 1;
        Reservoir(ires).CurrentAccPerRoute(i_r) = Reservoir(ires).CurrentAccPerRoute(i_r) + 1;
        NextEvent.ResList = [NextEvent.ResList ires];
        
        % Update entering curve
        Reservoir(ires).EntryTimes = [Reservoir(ires).EntryTimes CurrentTime];
        Reservoir(ires).EntryTimesPerRoute{i_r} = [Reservoir(ires).EntryTimesPerRoute{i_r} CurrentTime];
        Reservoir(ires).CurrentNinPerRoute(i_r) = Reservoir(ires).CurrentNinPerRoute(i_r) + 1;
        
        % Update reservoir mean speed
        Temp_nr = Reservoir(ires).CurrentAcc;
        Temp_param = Reservoir(ires).MFDfctParam;
        Temp_scf = Simulation.TripbasedSimuFactor;
        if Temp_nr == 0
            Temp_Vr = Reservoir(ires).FreeflowSpeed;
        else
            Temp_Vr = MFDfct(Temp_nr,Temp_scf*Temp_param)/Temp_nr;
        end
        Reservoir(ires).CurrentMeanSpeed = Temp_Vr;
        Reservoir(ires).MeanSpeed2(itime) = Temp_Vr;
        
        % Global waiting list for exiting the reservoir: order of exit
        Temp_vehlist = Reservoir(ires).VehList;
        if isempty(Temp_vehlist)
            Reservoir(ires).VehList = iveh;
        else
            i_p = Vehicle(iveh).PathIndex;
            Temp_dist = Vehicle(iveh).TripLength(i_p) - Vehicle(iveh).TraveledDistance(i_p); % remaining distance
            Temp_vehlist = [Temp_vehlist iveh]; % veh in last position by default
            Temp_vehlist2 = Temp_vehlist;
            i = length(Temp_vehlist2);
            if i > 1
                i_p = Vehicle(Temp_vehlist2(i-1)).PathIndex;
                Temp_dist2 = Vehicle(Temp_vehlist2(i-1)).TripLength(i_p) - Vehicle(Temp_vehlist2(i-1)).TraveledDistance(i_p); % remaining distance
                while i > 1 && Temp_dist < Temp_dist2
                    Temp_vehlist(i) = Temp_vehlist(i-1);
                    Temp_vehlist(i-1) = iveh;
                    i = i - 1;
                    if i > 1
                        i_p = Vehicle(Temp_vehlist2(i-1)).PathIndex;
                        Temp_dist2 = Vehicle(Temp_vehlist2(i-1)).TripLength(i_p) - Vehicle(Temp_vehlist2(i-1)).TraveledDistance(i_p); % remaining distance
                    end
                end
            end
            Reservoir(ires).VehList = Temp_vehlist;
        end
        
        % Set the first vehicle to exit for each route in r
        if Reservoir(ires).FirstVehPerRoute(i_r) == 0
            Reservoir(ires).FirstVehPerRoute(i_r) = iveh;
        end
        
    else % The vehicle has arrived at destination
        
        % Travel time of the vehicle
        Temp_tentry = Vehicle(iveh).EntryTimes(1);
        Temp_texit = Vehicle(iveh).ExitTimes(end);
        Route(iroute).TravelTime2(itime) = Temp_texit - Temp_tentry;
        
    end
    
    
    % Desired exit time in each reservoir
    %------------------------------------
    for r = NextEvent.ResList % loop only on reservoirs where accumulation changed
        Temp_nr = Reservoir(r).CurrentAcc;
        Temp_Pc = Reservoir(r).MaxProd;
        Temp_scf = Simulation.TripbasedSimuFactor;
        for i_r = Reservoir(r).ExitRoutesIndex % loop on all routes exiting Rr
            iveh = Reservoir(r).FirstVehPerRoute(i_r);
            if iveh > 0
                i_p = Vehicle(iveh).PathIndex;
                if Reservoir(r).CurrentAcc > Temp_scf*Reservoir(r).CritAcc && ...
                        (strcmp(Simulation.DivergeModel,'maxdem') || strcmp(Simulation.DivergeModel,'queuedyn'))
                    % Force the veh to exit in case of congested state
                    Vehicle(iveh).TraveledDistance(i_p) = Vehicle(iveh).TripLength(i_p);
                end
                if Reservoir(r).CurrentMeanSpeed > 0
                    Temp_time = CurrentTime + (Vehicle(iveh).TripLength(i_p) - Vehicle(iveh).TraveledDistance(i_p))/Reservoir(r).CurrentMeanSpeed;
                else
                    %Temp_time = Inf;
                    Temp_time = CurrentTime;
                end
                Temp_nrp = Reservoir(r).CurrentAccPerRoute(i_r);
                Temp_Lrp = Reservoir(r).TripLengthPerRoute(i_r);
                %Temp_time2 = Reservoir(r).LastExitTimePerRoute(i_r) + Temp_nr/Temp_nrp*Temp_Lrp/(Temp_scf*Temp_Pc);
                Temp_time2 = 0;
                Temp_time = max([Temp_time CurrentTime Temp_time2]);
                Reservoir(r).DesiredExitTimePerRoute(i_r) = Temp_time;
                Reservoir(r).DesiredExitVehPerRoute(i_r) = iveh;
            else
                Reservoir(r).DesiredExitTimePerRoute(i_r) = Inf;
                Reservoir(r).DesiredExitVehPerRoute(i_r) = 1;
            end
        end
        for i_r = Reservoir(r).DestRoutesIndex % loop on all routes ending in Rr
            iveh = Reservoir(r).FirstVehPerRoute(i_r);
            if iveh > 0
                i_p = Vehicle(iveh).PathIndex;
                if Reservoir(r).CurrentMeanSpeed > 0
                    Temp_time = CurrentTime + (Vehicle(iveh).TripLength(i_p) - Vehicle(iveh).TraveledDistance(i_p))/Reservoir(r).CurrentMeanSpeed;
                else
                    %Temp_time = Inf;
                    Temp_time = CurrentTime;
                end
                Temp_nrp = Reservoir(r).CurrentAccPerRoute(i_r);
                Temp_Lrp = Reservoir(r).TripLengthPerRoute(i_r);
                %Temp_time2 = Reservoir(r).LastExitTimePerRoute(i_r) + Temp_nr/Temp_nrp*Temp_Lrp/(Temp_scf*Temp_Pc);
                Temp_time2 = 0;
                Temp_time = max([Temp_time CurrentTime Temp_time2]);
                Reservoir(r).DesiredExitTimePerRoute(i_r) = Temp_time;
                Reservoir(r).DesiredExitVehPerRoute(i_r) = iveh;
            else
                Reservoir(r).DesiredExitTimePerRoute(i_r) = Inf;
                Reservoir(r).DesiredExitVehPerRoute(i_r) = 1;
            end
        end
    end
    
    
    % Desired entry time in each reservoir
    %-------------------------------------
    % Update the desired entry time if the current event is an entry (vehicle creation)
    if NextEvent.Type == 1
        iveh = NextEvent.VehID;
        ires = Vehicle(iveh).CurrentResID;
        iroute = Vehicle(iveh).RouteID;
        i_r = Route(iroute).ResRouteIndex(ires);
        
        Reservoir(ires).DemandTimeIndexPerRoute(i_r) = Reservoir(ires).DemandTimeIndexPerRoute(i_r) + 1;
        Temp_index = Reservoir(ires).DemandTimeIndexPerRoute(i_r);
        if Temp_index <= length(Reservoir(ires).DemandEntryTimePerRoute{i_r})
            Reservoir(ires).DesiredEntryTimePerRoute(i_r) = Reservoir(ires).DemandEntryTimePerRoute{i_r}(Temp_index);
            Reservoir(ires).DesiredEntryVehPerRoute(i_r) = Reservoir(ires).DemandEntryVehPerRoute{i_r}(Temp_index);
        else
            Reservoir(ires).DesiredEntryTimePerRoute(i_r) = Inf;
            Reservoir(ires).DesiredEntryVehPerRoute(i_r) = 1;
        end
        Reservoir(ires).DemandTimeIndex = Reservoir(ires).DemandTimeIndex + 1;
        Temp_index = Reservoir(ires).DemandTimeIndex;
        if Temp_index <= length(Reservoir(ires).DemandEntryTime)
            Reservoir(ires).DesiredEntryTime = Reservoir(ires).DemandEntryTime(Temp_index); % only for route beginning
            Reservoir(ires).DesiredEntryVeh = Reservoir(ires).DemandEntryVeh(Temp_index); % only for route beginning
        else
            Reservoir(ires).DesiredEntryTime = Inf;
            Reservoir(ires).DesiredEntryVeh = 1;
        end
    end
    % Update the desired entry time with the desired exit time from previous res (route transfer)
    Temp_reslist = [];
    for r = NextEvent.ResList % loop only on reservoirs where accumulation changed
        for i_r = Reservoir(r).ExitRoutesIndex % loop on all routes exiting Rr
            iroute = Reservoir(r).RoutesID(i_r);
            if r ~= Route(iroute).ResDestinationID
                r2 = Route(iroute).ResPath(Reservoir(r).RoutesPathIndex(i_r)+1); % next reservoir in the route
                i_r2 = Route(iroute).ResRouteIndex(r2);
                Reservoir(r2).DesiredEntryTimePerRoute(i_r2) = Reservoir(r).DesiredExitTimePerRoute(i_r);
                Reservoir(r2).DesiredEntryVehPerRoute(i_r2) = Reservoir(r).DesiredExitVehPerRoute(i_r);
                Temp_reslist = [Temp_reslist r2];
            end
        end
    end
    
    
    % Entry supply time per route in each reservoir
    %----------------------------------------------
    for r = NextEvent.ResList % loop only on reservoirs where accumulation changed
        
        Temp_Nroutes = length(Reservoir(r).RoutesID);
        Temp_Nsmooth = sum(Reservoir(r).CurrentAccPerRoute > 0);
        %Temp_Nsmooth = length(Reservoir(r).RoutesID) + 2;
        Temp_qinr2 = zeros(1,Temp_Nroutes);
        for i_r = 1:Temp_Nroutes
            if Reservoir(r).DesiredEntryTimePerRoute(i_r) > Reservoir(r).LastEntryTimePerRoute(i_r) ...
                    && Reservoir(r).LastEntryTimePerRoute(i_r) > -Inf && Reservoir(r).DesiredEntryTimePerRoute(i_r) < Inf
                Temp_qinr2(i_r) = 1/(Reservoir(r).DesiredEntryTimePerRoute(i_r) - Reservoir(r).LastEntryTimePerRoute(i_r));
            elseif Reservoir(r).DesiredEntryTimePerRoute(i_r) == Inf
                Temp_qinr2(i_r) = 0;
            else
                iroute = Reservoir(r).RoutesID(i_r);
                inode = Route(iroute).NodePath(Reservoir(r).RoutesPathIndex(i_r)); % entry/origin node for iroute in Rr
                Temp_qinr2(i_r) = Simulation.TripbasedSimuFactor.*MacroNode(inode).Supply;
            end
        end
        
        % Effective inflow for internal origins (not restricted)
        Reservoir(r).InternalProd = 0;
        for i_n = Reservoir(r).OriNodesIndex % loop on all origin nodes in Rr
            inode = Reservoir(r).MacroNodesID(i_n);
            Temp_indexes = Reservoir(r).NodeRoutesIndex{i_n};
            if ~isempty(Temp_indexes)
                Temp_dem = Temp_qinr2(Temp_indexes);
                Temp_demtot = sum(Temp_dem);
                if Temp_demtot > 0
                    Temp_mergecoeff = Temp_dem./Temp_demtot;
                else
                    Temp_mergecoeff = ones(1,length(Temp_indexes));
                end
                Temp_sup = Simulation.TripbasedSimuFactor.*MacroNode(inode).Supply;
                Temp_inflow = mergeFair(Temp_dem,Temp_sup,Temp_mergecoeff);
                Temp_supplytimes = Reservoir(r).LastEntryTimePerRoute(Temp_indexes) + 1./Temp_inflow;
                Reservoir(r).EntrySupplyTimePerRoute(Temp_indexes) = max([Temp_supplytimes; CurrentTime*ones(1,length(Temp_indexes))]);
                Temp_Ltrip = Reservoir(r).TripLengthPerRoute(Temp_indexes);
                Reservoir(r).InternalProd = Reservoir(r).InternalProd + sum(Temp_Ltrip.*Temp_inflow);
            end
        end
        
        % Entry merge coefficients for entering routes
        if strcmp(Simulation.MergeModel,'endogenous')
            % Endogenous merge (for entering productions)
            Temp_indexes = Reservoir(r).EntryRoutesIndex;
            Temp_nrp = Reservoir(r).CurrentAccPerRoute(Temp_indexes);
            Temp_nr_entry = sum(Temp_nrp);
            if Temp_nr_entry > 0
                Reservoir(r).MergeCoeffPerRoute(Temp_indexes) = (Temp_nrp > 0).*Temp_nrp./Temp_nr_entry + (Temp_nrp <= 0).*1;
            end
        elseif strcmp(Simulation.MergeModel,'demprorata') || strcmp(Simulation.MergeModel,'demfifo')
            % Demand pro-rata flow merge
            Temp_indexes = Reservoir(r).EntryRoutesIndex;
            Temp_dem = Temp_qinr2(Temp_indexes);
            Temp_demtot = sum(Temp_dem);
            if Temp_demtot > 0
                Reservoir(r).MergeCoeffPerRoute(Temp_indexes) = Temp_dem./Temp_demtot;
            end
        elseif strcmp(Simulation.MergeModel,'equiproba')
            % Equi-probability for all transfer inflows
            Temp_indexes = Reservoir(r).EntryRoutesIndex;
            Temp_Nroutes = length(Temp_indexes);
            if Temp_Nroutes > 0
                Reservoir(r).MergeCoeffPerRoute(Temp_indexes) = ones(1,Temp_Nroutes)./Temp_Nroutes;
            end
        end
        
        % Inflow limitation due to node supply at entry (border supply)
        for i_n = Reservoir(r).EntryNodesIndex % loop on all entry border nodes in Rr
            inode = Reservoir(r).MacroNodesID(i_n);
            Temp_indexes = Reservoir(r).NodeRoutesIndex{i_n};
            if ~isempty(Temp_indexes)
                Temp_sup = Simulation.TripbasedSimuFactor.*MacroNode(inode).Supply;
                % Modify the original desired entry times for the endogenous merge
                if strcmp(Simulation.MergeModel,'endogenous')
                    Temp_demtimes = Reservoir(r).DesiredEntryTimePerRoute(Temp_indexes);
                    Temp_lasttimes = Reservoir(r).LastEntryTimePerRoute(Temp_indexes);
                    Temp_mergecoeff = Reservoir(r).MergeCoeffPerRoute(Temp_indexes);
                    Temp_mergecoefftot = sum(Temp_mergecoeff);
                    Temp_entrytimes = mergetimeFair(ones(1,length(Temp_indexes)),Temp_demtimes,Temp_lasttimes,1,Temp_sup,Temp_mergecoeff./Temp_mergecoefftot);
                    Reservoir(r).DesiredEntryTimePerRoute(Temp_indexes) = max([Temp_entrytimes; Reservoir(r).DesiredEntryTimePerRoute(Temp_indexes)]);
                end
                % Use directly the node supply for other merge models
                Reservoir(r).EntryNodeSupplyPerRoute(Temp_indexes) = Temp_sup;
            end
        end
        
        % Effective inflow for entering routes
        if strcmp(Simulation.MergeModel,'endogenous')
            % Endogenous merge (for entering productions)
            Temp_nr = Reservoir(r).CurrentAcc;
            Temp_scf = Simulation.TripbasedSimuFactor;
            Temp_param = Reservoir(r).EntryfctParam;
            Temp_prodsupply = Entryfct(Temp_nr,Temp_scf*Temp_param) - Reservoir(r).InternalProd; % SCF already included in InternalProd
            Temp_indexes = Reservoir(r).EntryRoutesIndex;
            Temp_Lrp = Reservoir(r).TripLengthPerRoute(Temp_indexes);
            Temp_mergecoeff = Reservoir(r).MergeCoeffPerRoute(Temp_indexes);
            Temp_demtimes = Reservoir(r).DesiredEntryTimePerRoute(Temp_indexes);
            Temp_lasttimes = Reservoir(r).LastEntryTimePerRoute(Temp_indexes);
            
            Temp_supplytimes = mergetimeFair(Temp_Lrp,Temp_demtimes,Temp_lasttimes,1,Temp_prodsupply,Temp_mergecoeff);
            Reservoir(r).EntrySupplyTimePerRoute(Temp_indexes) = max([Temp_supplytimes; CurrentTime*ones(1,length(Temp_indexes))]);
        else
            % Other merge models (for inflows)
            Temp_nr = Reservoir(r).CurrentAcc;
            Temp_scf = Simulation.TripbasedSimuFactor;
            Temp_param = Reservoir(r).EntryfctParam;
            Temp_prodsupply = Entryfct(Temp_nr,Temp_scf*Temp_param) - Reservoir(r).InternalProd; % SCF already included in InternalProd
            Temp_indexes = Reservoir(r).EntryRoutesIndex;
            Temp_Lrp = Reservoir(r).TripLengthPerRoute(Temp_indexes);
            Temp_proddem = Temp_Lrp.*Temp_qinr2(Temp_indexes);
            Temp_mergecoeff = Reservoir(r).MergeCoeffPerRoute(Temp_indexes);
            Temp_nrp = Reservoir(r).CurrentAccPerRoute(Temp_indexes);
            Temp_nr_entry = sum(Temp_nrp);
            if Temp_nr_entry > 0
                Temp_Lr_entry = Temp_nr_entry/sum(Temp_nrp./Temp_Lrp);
            else
                Temp_Lr_entry = mean(Temp_Lrp);
            end
            if sum(Temp_proddem) < Temp_prodsupply
                Temp_flowsupply = Inf;
            else
                Temp_flowsupply = Temp_prodsupply/Temp_Lr_entry;
            end
            Temp_demtimes = Reservoir(r).DesiredEntryTimePerRoute(Temp_indexes);
            Temp_lasttimes = Reservoir(r).LastEntryTimePerRoute(Temp_indexes);
            
            Temp_supplytimes = mergetimeFair(ones(1,length(Temp_indexes)),Temp_demtimes,Temp_lasttimes,1,Temp_flowsupply,Temp_mergecoeff);
            Reservoir(r).EntrySupplyTimePerRoute(Temp_indexes) = max([Temp_supplytimes; CurrentTime*ones(1,length(Temp_indexes))]);
            Reservoir(r).EntrySupplyTime = Reservoir(r).LastEntryTime + 1/Temp_flowsupply; % global entry supply time
        end
    end
    
    
    % Reservoir exit supply limitation (destination)
    %---------------------------------
    for r = Global.DestResList % loop on all destination reservoirs
        
        % Outflow demand estimation
        Temp_scf = Simulation.TripbasedSimuFactor;
        Temp_Nroutes = length(Reservoir(r).RoutesID);
        Temp_qoutr = zeros(1,Temp_Nroutes);
        for i_r = 1:Temp_Nroutes
            if Reservoir(r).DesiredExitTimePerRoute(i_r) > Reservoir(r).LastExitTimePerRoute(i_r) ...
                    && Reservoir(r).LastExitTimePerRoute(i_r) > -Inf && Reservoir(r).DesiredExitTimePerRoute(i_r) < Inf
                Temp_qoutr(i_r) = 1/(Reservoir(r).DesiredExitTimePerRoute(i_r) - Reservoir(r).LastExitTimePerRoute(i_r));
            elseif Reservoir(r).DesiredExitTimePerRoute(i_r) == Inf
                Temp_qoutr(i_r) = 0;
            else
                Temp_nr = Reservoir(r).CurrentAcc;
                Temp_nrp = Reservoir(r).CurrentAccPerRoute(i_r);
                Temp_Lrp = Reservoir(r).TripLengthPerRoute(i_r);
                Temp_Pc = Reservoir(r).MaxProd;
                Temp_qoutr(i_r) = Temp_nrp/Temp_nr*(Temp_scf*Temp_Pc)/Temp_Lrp; % maximum outflow
            end
        end
        
        % Exit merge coefficients (merge for several outflows ending in the same node)
        Temp_dem = Temp_qoutr;
        Temp_demtot = sum(Temp_dem);
        if Temp_demtot > 0
            Reservoir(r).ExitCoeffPerRoute = Temp_dem./Temp_demtot;
        end
        
        % Exit supply for internal destinations
        for i_n = Reservoir(r).DestNodesIndex % loop on all destination nodes in Rr
            inode = Reservoir(r).MacroNodesID(i_n);
            Temp_indexes = Reservoir(r).NodeRoutesIndex{i_n};
            if ~isempty(Temp_indexes)
                Temp_mergecoeff = Reservoir(r).ExitCoeffPerRoute(Temp_indexes);
                Temp_mergecoefftot = sum(Temp_mergecoeff);
                Temp_demtimes = Reservoir(r).DesiredExitTimePerRoute(Temp_indexes);
                Temp_lasttimes = Reservoir(r).LastExitTimePerRoute(Temp_indexes);
                Temp_sup = Simulation.TripbasedSimuFactor.*MacroNode(inode).Supply;
                
                Temp_supplytimes = mergetimeFair(ones(1,length(Temp_indexes)),Temp_demtimes,Temp_lasttimes,1,Temp_sup,Temp_mergecoeff./Temp_mergecoefftot);
                Reservoir(r).ExitSupplyTimePerRoute(Temp_indexes) = max([Temp_supplytimes; CurrentTime*ones(1,length(Temp_indexes))]);
            end
        end
    end
    
    
    % Reservoir exit supply time (in the middle of a route)
    % Next entry in reservoirs
    %---------------------------
    NextEvent.ExitResList = Global.DestResList;
    for r = NextEvent.ResList % loop only on reservoirs where accumulation changed
        
        % Exit supply for external destinations and transfer to another reservoir
        for i_n = Reservoir(r).ExitNodesIndex % loop on all exit border nodes in Rr
            inode = Reservoir(r).MacroNodesID(i_n);
            Temp_indexes = Reservoir(r).NodeRoutesIndex{i_n};
            if ~isempty(Temp_indexes)
                Temp_mergecoeff = Reservoir(r).ExitCoeffPerRoute(Temp_indexes);
                Temp_mergecoefftot = sum(Temp_mergecoeff);
                Temp_demtimes = Reservoir(r).DesiredExitTimePerRoute(Temp_indexes);
                Temp_lasttimes = Reservoir(r).LastExitTimePerRoute(Temp_indexes);
                Temp_sup = Simulation.TripbasedSimuFactor.*MacroNode(inode).Supply;
                
                Temp_supplytimes = mergetimeFair(ones(1,length(Temp_indexes)),Temp_demtimes,Temp_lasttimes,1,Temp_sup,Temp_mergecoeff./Temp_mergecoefftot);
                i = 1;
                for i_r = Temp_indexes
                    iroute = Reservoir(r).RoutesID(i_r);
                    if r == Route(iroute).ResDestinationID
                        Reservoir(r).ExitSupplyTimePerRoute(i_r) = max([Temp_supplytimes(i) CurrentTime]);
                    else
                        r2 = Route(iroute).ResPath(Reservoir(r).RoutesPathIndex(i_r)+1); % next reservoir in the route
                        i_r2 = Route(iroute).ResRouteIndex(r2);
                        if strcmp(Simulation.MergeModel,'endogenous')
                            Reservoir(r).ExitSupplyTimePerRoute(i_r) = Reservoir(r2).EntrySupplyTimePerRoute(i_r2); % node supply already included in the inflow supply
                        else
                            Temp_time = Reservoir(r2).LastEntryTime + 1/Reservoir(r2).EntryNodeSupplyPerRoute(i_r2);
                            Reservoir(r).ExitSupplyTimePerRoute(i_r) = max([Reservoir(r2).EntrySupplyTime Temp_time]);
                        end
                        NextEvent.ExitResList = [NextEvent.ExitResList r];
                    end
                    i = i + 1;
                end
            end
        end
        
        % Next entry in the current reservoir (beginning of a route)
        Temp_indexes = [];
        for i_r = 1:length(Reservoir(r).RoutesID)
            iroute = Reservoir(r).RoutesID(i_r);
            if r == Route(iroute).ResOriginID
                Temp_indexes = [Temp_indexes i_r];
            end
        end
        if ~isempty(Temp_indexes)
            if strcmp(Simulation.MergeModel,'demfifo') && length(Reservoir(r).EntryNodesIndex) == 1
                % FIFO discipline merge only if one external entry node
                Temp_demandtime = Reservoir(r).DesiredEntryTime; % not necessarily route beginning!
                Temp_supplytime = Reservoir(r).EntrySupplyTime;
                iveh = Reservoir(r).DesiredEntryVeh;
                iroute = Vehicle(iveh).RouteID;
                i_r = Route(iroute).ResRouteIndex(r);
                Temp_time2 = Reservoir(r).LastEntryTime + 1/Reservoir(r).EntryNodeSupplyPerRoute(i_r);
                Reservoir(r).NextEntryTime = max([Temp_demandtime Temp_supplytime Temp_time2]);
                Reservoir(r).NextEntryVehID = iveh;
            else
                % Other merge models
                Temp_demandtimes = Reservoir(r).DesiredEntryTimePerRoute(Temp_indexes);
                Temp_supplytimes = Reservoir(r).EntrySupplyTimePerRoute(Temp_indexes);
                [Temp_time, i_r] = selectmintime(Temp_demandtimes,Temp_supplytimes);
                iveh = Reservoir(r).DesiredEntryVehPerRoute(Temp_indexes(i_r));
                if strcmp(Simulation.MergeModel,'endogenous')
                    Temp_time2 = 0;
                else
                    Temp_time2 = Reservoir(r).LastEntryTime + 1/Reservoir(r).EntryNodeSupplyPerRoute(Temp_indexes(i_r));
                end
                Reservoir(r).NextEntryTime = max([Temp_time Temp_time2 CurrentTime]);
                Reservoir(r).NextEntryVehID = iveh;
            end
        else
            Reservoir(r).NextEntryTime = Inf;
            Reservoir(r).NextEntryVehID = 1;
        end
    end
    NextEvent.ExitResList = unique(NextEvent.ExitResList);
    
    
    % Next exit in reservoirs
    %------------------------
    Temp_reslist = unique([NextEvent.ResList NextEvent.ExitResList]);
    for r = Temp_reslist
        % Next exit from the current reservoir
        if strcmp(Simulation.DivergeModel,'maxdem') || strcmp(Simulation.DivergeModel,'queuedyn')
            if ~isempty(Reservoir(r).VehList)
                iveh = Reservoir(r).VehList(1);
                iroute = Vehicle(iveh).RouteID;
                i_r = Route(iroute).ResRouteIndex(r);
                Temp_nr = Reservoir(r).CurrentAcc;
                Temp_Pc = Reservoir(r).MaxProd;
                Temp_Lrp = Reservoir(r).TripLengthPerRoute;
                Temp_nrp = Reservoir(r).CurrentAccPerRoute;
                if Temp_nr == 0
                    Temp_Lr = mean(Temp_Lrp);
                else
                    Temp_Lr = Temp_nr/sum(Temp_nrp./Temp_Lrp);
                end
                Temp_time2 = Reservoir(r).LastExitTime + Temp_Lr/(Temp_scf*Temp_Pc);
                Temp_demandtime = Reservoir(r).DesiredExitTimePerRoute(i_r);
                Temp_supplytime = Reservoir(r).ExitSupplyTimePerRoute(i_r);
                Temp_time = max([Temp_demandtime Temp_supplytime Temp_time2 CurrentTime]);
            else
                iveh = 0;
                Temp_time = Inf;
            end
            Reservoir(r).NextExitTime = Temp_time;
            Reservoir(r).NextExitVehID = iveh;
        else
            Temp_demandtimes = Reservoir(r).DesiredExitTimePerRoute;
            Temp_supplytimes = Reservoir(r).ExitSupplyTimePerRoute;
            [Temp_time, Temp_iroute] = selectmintime(Temp_demandtimes,Temp_supplytimes);
            Reservoir(r).NextExitTime = max([Temp_time CurrentTime]);
            Reservoir(r).NextExitVehID = Reservoir(r).DesiredExitVehPerRoute(Temp_iroute);
        end
    end
    
    
    % Next global entry and exit
    %---------------------------
    Temp_entrytimes = Inf*ones(1,NumRes);
    Temp_exittimes = Inf*ones(1,NumRes);
    for r = 1:NumRes
        if ~isempty(Reservoir(r).RoutesID)
            Temp_entrytimes(r) = Reservoir(r).NextEntryTime;
            Temp_exittimes(r) = Reservoir(r).NextExitTime;
        end
    end
    indexsetmin = find(Temp_entrytimes == min(Temp_entrytimes));
    Nset = length(indexsetmin);
    % Uniform random draw of an entry among the min times
    i = randi(Nset);
    ires = indexsetmin(i);
    NextEntry.Time = Reservoir(ires).NextEntryTime;
    NextEntry.VehID = Reservoir(ires).NextEntryVehID;
    
    indexsetmin = find(Temp_exittimes == min(Temp_exittimes));
    Nset = length(indexsetmin);
    % Uniform random draw of an entry among the min times
    i = randi(Nset);
    ires = indexsetmin(i);
    NextExit.Time = Reservoir(ires).NextExitTime;
    NextExit.VehID = Reservoir(ires).NextExitVehID;
    
    
    % Next event
    %-----------
    NextEvent.Time = min([NextEntry.Time NextExit.Time]);
    NextEvent.ElapsedTime = NextEvent.Time - CurrentTime;
    
    if NextEvent.Time == NextEntry.Time
        % If the next event is an entry
        iveh = NextEntry.VehID;
        ires = Vehicle(iveh).CurrentResID;
        
        NextEvent.VehID = iveh;
        NextEvent.Type = 1;
        
        % Next entry time for the current reservoir
        Reservoir(ires).NextEntryTime = Inf;
        
    else
        % If the next event is an exit
        iveh = NextExit.VehID;
        ires = Vehicle(iveh).CurrentResID;
        
        NextEvent.VehID = iveh;
        NextEvent.Type = 2;
        
        % Next exit time for the current reservoir
        Reservoir(ires).NextExitTime = Inf;
        
    end
    
%     r = 2;
%     if 2500 <= CurrentTime && CurrentTime < 4500
%         fprintf('%s%3.1f \t %s%3.1f \t %s%i \t %s%i \t %s%i \n','time=',CurrentTime,'nexteventtime=',NextEvent.Time,'r=',Vehicle(NextEvent.VehID).CurrentResID,'iveh=',NextEvent.VehID,'route=',Vehicle(NextEvent.VehID).RouteID)
%         fprintf('%s \n','Last entry time per route:')
%         for i_r = 1:length(Reservoir(r).RoutesID)
%             fprintf('%3.1f \t',Reservoir(r).LastEntryTimePerRoute(i_r))
%         end
%         fprintf('\n%s \n','Entry supply time per route:')
%         for i_r = 1:length(Reservoir(r).RoutesID)
%             fprintf('%3.1f \t',Reservoir(r).EntrySupplyTimePerRoute(i_r))
%         end
%         fprintf('\n%s \n','Desired exit time per route:')
%         for i_r = 1:length(Reservoir(r).RoutesID)
%             fprintf('%3.1f \t',Reservoir(r).DesiredExitTimePerRoute(i_r))
%         end
%         fprintf('\n%s \n','Exit supply time per route:')
%         for i_r = 1:length(Reservoir(r).RoutesID)
%             fprintf('%3.1f \t',Reservoir(r).ExitSupplyTimePerRoute(i_r))
%         end
%         fprintf('\n%s%3.1f \t %s%3.1f \n','nextentry=',Reservoir(r).NextEntryTime,'nextexit=',Reservoir(r).NextExitTime)
%         fprintf('%s\n',' ')
%     end
    
    %keyboard
    
    % Time update
    %------------
    ElapsedTime = NextEvent.ElapsedTime;
    CurrentTime = NextEvent.Time;
    itime = itime + 1;
    
end



%% Post-processing
%--------------------------------------------------------------------------

Global.NumEvents = itime - 1;

Global.SimulTime = Global.SimulTime(1:Global.NumEvents);
Global.VehID = Global.VehID(1:Global.NumEvents);

% Travel time resampling
for iroute = 1:NumRoutes
    Route(iroute).TravelTime2 = Route(iroute).TravelTime2(1:Global.NumEvents);
end
for iroute = 1:NumRoutes
    for i = (Global.CurrentTimeID+1):Global.NumEvents
        if Route(iroute).TravelTime2(i) == 0
            Route(iroute).TravelTime2(i) = Route(iroute).TravelTime2(i-1);
        end
    end
    Temp_t0 = Global.SimulTime(1:Global.NumEvents);
    Temp_TT0 = Route(iroute).TravelTime2(1:Global.NumEvents);
    Temp_t = Simulation.Time(1:Temp_EndTimeID);
    Temp_TT = resamp(Temp_t,Temp_t0,Temp_TT0);
    Route(iroute).TravelTime(Temp_StartTimeID:Temp_EndTimeID) = Temp_TT(Temp_StartTimeID:Temp_EndTimeID);
    for i = Temp_StartTimeID:Temp_EndTimeID
        if Route(iroute).TravelTime(i) == 0
            Route(iroute).TravelTime(i) = Route(iroute).FreeFlowTravelTime;
        end
    end
    if Temp_StartTimeID > 1
        % to ensure continuity on the evolution
        Route(iroute).TravelTime(Temp_StartTimeID) = Route(iroute).TravelTime(Temp_StartTimeID-1);
    end
    Route(iroute).TravelTime(NumTimes) = Route(iroute).TravelTime(NumTimes-1);
end

% Mean speed resampling
for r = 1:NumRes
    Reservoir(r).MeanSpeed2 = Reservoir(r).MeanSpeed2(1:Global.NumEvents);
end
for r = 1:NumRes
    for i = (Global.CurrentTimeID+1):Global.NumEvents
        if Reservoir(r).MeanSpeed2(i) == 0
            Reservoir(r).MeanSpeed2(i) = Reservoir(r).MeanSpeed2(i-1);
        end
    end
    Temp_t0 = Global.SimulTime(1:Global.NumEvents);
    Temp_V0 = Reservoir(r).MeanSpeed2(1:Global.NumEvents);
    Temp_t = Simulation.Time(1:Temp_EndTimeID);
    Temp_V = resamp(Temp_t,Temp_t0,Temp_V0);
    Reservoir(r).MeanSpeed(Temp_StartTimeID:Temp_EndTimeID) = Temp_V(Temp_StartTimeID:Temp_EndTimeID);
    if Temp_StartTimeID > 1
        % to ensure continuity on the evolution
        Reservoir(r).MeanSpeed(Temp_StartTimeID) = Reservoir(r).MeanSpeed(Temp_StartTimeID-1);
    end
end





