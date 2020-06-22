%% ROUTE CALCULATION
%--------------------------------------------------------------------------

% Simulation time
Simulation.Time = 0:Simulation.TimeStep:Simulation.Duration; % vector of times [s]
SimulationDuration = Simulation.Duration;
TimeStep = Simulation.TimeStep;
CurrentTime = Simulation.Time;
NumTimes = floor(SimulationDuration/TimeStep) + 1; % Number of times
NumPaths = Assignment.NumShortestPath;

% Set default parameters (only those not set in SimulSettings already)
SetDefaultParam

% Calculate the number of modes
Temp_carmode = 0;
Temp_PTmode = 0;


if Assignment.PredefRoute == 0
    % If the routes are calculated by the assignment procedure
    %---------------------------------------------------------
    
    iroute = 1;
    
    for od = 1:NumODmacro % loop on all OD
        
        o = ODmacro(od).ResOriginID;
        d = ODmacro(od).ResDestinationID;
        ODmacro(od).RoutesID = [];
        
        Temp_demdata = 0;
        for k = 1:length(ODmacro(od).Demand)
            Temp_demdata = Temp_demdata + sum(ODmacro(od).Demand(k).Data);
        end
        
        if Temp_demdata > 0 % if the OD (o,d) is used at one time at least
            
            % Calculate cost for each possible route of OD (o,d)
            Temp_TTroutes = zeros(1,ODmacro(od).NumPossibleRoutes);
            for iroute0 = 1:ODmacro(od).NumPossibleRoutes % loop on all possible routes for OD (o,d)
                i_r = 1;
                for r = ODmacro(od).PossibleRoute(iroute0).ResPath % loop on all reservoirs in possible route iroute0
                    % Add cost of reservoir r
                    if Assignment.CurrentPeriodID == 1
                        Temp_TTroutes(iroute0) = Temp_TTroutes(iroute0) + ...
                            ODmacro(od).PossibleRoute(iroute0).TripLengths(i_r)/Reservoir(r).FreeflowSpeed;
                    else
                        Temp_TTroutes(iroute0) = Temp_TTroutes(iroute0) + ...
                            ODmacro(od).PossibleRoute(iroute0).TripLengths(i_r)/mean(Reservoir(r).MeanSpeed);
                    end
                    i_r = i_r + 1;
                end
            end
            
            % Sort all the possible routes by ascending cost
            [Temp_TTroutes, Temp_shortestpath] = sort(Temp_TTroutes);
            
            % Save the routes only for the desired number of shortest paths
            Temp_paths = min([NumPaths ODmacro(od).NumPossibleRoutes]);
            for p = 1:Temp_paths
                Route(iroute).ODmacroID = od;
                Route(iroute).ResPath = ODmacro(od).PossibleRoute(Temp_shortestpath(p)).ResPath;
                Route(iroute).NodePath = ODmacro(od).PossibleRoute(Temp_shortestpath(p)).NodePath;
                Route(iroute).TripLengths = ODmacro(od).PossibleRoute(Temp_shortestpath(p)).TripLengths;
                Route(iroute).NumMicroTrips = ODmacro(od).PossibleRoute(Temp_shortestpath(p)).NumMicroTrips;
                Route(iroute).ResOriginID = o;
                Route(iroute).ResDestinationID = d;
                Route(iroute).ResPathID = p;
                Route(iroute).NodeOriginID = Route(iroute).NodePath(1);
                Route(iroute).NodeDestinationID = Route(iroute).NodePath(end);
                Route(iroute).Length = sum(Route(iroute).TripLengths);
                Route(iroute).TotalTime = Temp_TTroutes(p);
                Route(iroute).FreeFlowTravelTime = Temp_TTroutes(p);
                Route(iroute).OldTT = Route(iroute).TotalTime;
                
                ODmacro(od).RoutesID = [ODmacro(od).RoutesID iroute];
                
                iroute = iroute + 1;
            end
        end
    end
    
    % Number of routes
    NumRoutes = iroute - 1;
    
else
    % If the routes with their demand are given prior to the simulation
    %------------------------------------------------------------------
    for od = 1:NumODmacro % loop on all OD
        ODmacro(od).RoutesID = [];
    end
    
    NumRoutes = length(Route);
    
    for iroute = 1:NumRoutes
        
        od = Route(iroute).ODmacroID;
        o = Route(iroute).ResPath(1);
        d = Route(iroute).ResPath(end);
        p = length(ODmacro(od).RoutesID) + 1;
        Route(iroute).Mode = Route(iroute).Demand0(1).Purpose;
        if strcmp(Route(iroute).Demand0(j).Purpose,'cartrip')
            Route(iroute).ModeID = 1;
        elseif strcmp(Route(iroute).Demand0(j).Purpose,'publictransport')
            Route(iroute).ModeID = 2;
        end
        Route(iroute).ResOriginID = o;
        Route(iroute).ResDestinationID = d;
        Route(iroute).ResPathID = p;
        Route(iroute).NodeOriginID = Route(iroute).NodePath(1);
        Route(iroute).NodeDestinationID = Route(iroute).NodePath(end);
        Route(iroute).Length = sum(Route(iroute).TripLengths);
        Temp_TT = zeros(1,length(Route(iroute).TripLengths));
        i_p = 1;
        for r = Route(iroute).ResPath % loop on all reservoirs in path p
            ModeID = Route(iroute).ModeID;
            % Temp_TT = Temp_TT + Route(iroute).TripLengths(i_p)/Reservoir(r).FreeflowSpeed(ModeID);
            Temp_TT(i_p) = Route(iroute).TripLengths(i_p)/Reservoir(r).FreeflowSpeed(ModeID);
            i_p = i_p + 1;
        end
        Route(iroute).TotalTime = sum(Temp_TT); % route free-flow TT
        Route(iroute).FreeFlowTTPerReservoir = Temp_TT;
        Route(iroute).FreeFlowTravelTime = sum(Temp_TT);
        Route(iroute).OldTT = Route(iroute).TotalTime;
        
        ODmacro(od).RoutesID = [ODmacro(od).RoutesID iroute];
        
        % Append the route demand to the corresponding OD demand
        for j = 1:length(Route(iroute).Demand0)
            if strcmp(Route(iroute).Demand0(j).Purpose,'cartrip')
                Temp_carmode = 1;
                Temp_times = Route(iroute).Demand0(j).Time;
                Temp_data = Route(iroute).Demand0(j).Data;
                [Temp_times, Temp_data] = stairfct(Temp_times,Temp_data,TimeStep,0,SimulationDuration);
                
                i = 1;
                Route(iroute).Demand = Temp_data; % for the acc-based solver
                ODmacro(od).Demand(i).Time = Temp_times; % by default the time is discretized with the simul timestep
                ODmacro(od).Demand(i).Data = ODmacro(od).Demand(i).Data + Temp_data;
            end
            if strcmp(Route(iroute).Demand0(j).Purpose,'publictransport')
                Temp_PTmode = 1;
                Temp_times = Route(iroute).Demand0(j).Time;
                Temp_data = Route(iroute).Demand0(j).Data;
                [Temp_times, Temp_data] = stairfct(Temp_times,Temp_data,TimeStep,0,SimulationDuration);
                
                i = 2;
                Route(iroute).Demand = Temp_data; % for the acc-based solver
                ODmacro(od).Demand(i).Time = Temp_times; % by default the time is discretized with the simul timestep
                ODmacro(od).Demand(i).Data = ODmacro(od).Demand(i).Data + Temp_data;
            end
        end
    end
    Route = rmfield(Route,'Demand0'); % delete this field to save memory
end

% Number of modes
NumModes = Temp_carmode + Temp_PTmode;
Simulation.NumModes = NumModes;

% Calculating the total number of vehicules per OD, if the user wants to
% update the trip lengths for method 1 and 2.
for od = 1:NumODmacro
    if ~isempty(ODmacro(od).RoutesID)
        Temp_time = [ODmacro(od).Demand(1).Time SimulationDuration];
        Temp_demand = ODmacro(od).Demand(1).Data;
        
        Temp_nveh = 0;
        for k = 1:length(Temp_demand)
            Temp_nveh = Temp_nveh + Temp_demand(k)*(Temp_time(k+1) - Temp_time(k));
        end
        ODmacro(od).NVehicules = Temp_nveh;
    end
end

% Change trip lengths for Method 1 (one mean value per reservoir)
if Assignment.TripLengthMethod == 1
    Temp_NbTripMin = 1; % threshold to consider a macro-path or not (nb min of micro-trips)
    Temp_Ltrip = cell(1,NumRes);
    for od = 1:NumODmacro
        if ODmacro(od).NumPossibleRoutes > 0
            for iroute = 1:ODmacro(od).NumPossibleRoutes
                if ODmacro(od).PossibleRoute(iroute).NumMicroTrips > Temp_NbTripMin
                    i_r = 1;
                    for r = ODmacro(od).PossibleRoute(iroute).ResPath
                        Temp_Ltrip{r} = [Temp_Ltrip{r} ODmacro(od).PossibleRoute(iroute).TripLengths(i_r)];
                        i_r = i_r + 1;
                    end
                end
            end
        end
    end
    for iroute = 1:NumRoutes
        i_r = 1;
        for r = Route(iroute).ResPath
            Route(iroute).TripLengths(i_r) = mean(Temp_Ltrip{r});
            i_r = i_r + 1;
        end
    end
end


% Append to the Reservoir structure
for r = 1:NumRes
    Reservoir(r).TripLengthPerRoute = []; % trip lengths of the routes crossing the reservoir
    Reservoir(r).FreeFlowTTPerRoute = []; % free flow travel times of the routes crossing the reservoir
    Reservoir(r).RoutesID = []; % routes ID crossing the reservoir
    Reservoir(r).RoutesPathIndex = []; % indexes of the reservoir in the paths
    Reservoir(r).OriginRes = []; % reservoir origins of the routes
    Reservoir(r).DestinationRes = []; % reservoir destinations of the routes
    Reservoir(r).ODmacroID = cell(1,NumModes); % ODmacro IDs of the routes
    
    Reservoir(r).RoutesNodeID = []; % entry and exit nodes ID for each route crossing the reservoir
    Reservoir(r).NodeRoutesIndex = cell(NumModes,length(Reservoir(r).MacroNodesID)); % routes crossing each node in the reservoir, indexes for Reservoir(i).RoutesID
    
    Reservoir(r).OriRoutesIndex = cell(1,NumModes); % routes originating in the reservoir, indexes for Reservoir(i).RoutesID
    Reservoir(r).DestRoutesIndex = cell(1,NumModes); % routes ending in the reservoir, indexes for Reservoir(i).RoutesID
    Reservoir(r).EntryRoutesIndex = cell(1,NumModes); % routes entering the reservoir, indexes for Reservoir(i).RoutesID
    Reservoir(r).ExitRoutesIndex = cell(1,NumModes); % routes exiting the reservoir, indexes for Reservoir(i).RoutesID
    
    Reservoir(r).OriNodesIndex = cell(1,NumModes); % origin nodes in the reservoir, indexes for Reservoir(i).MacroNodesID
    Reservoir(r).DestNodesIndex = cell(1,NumModes); % destination nodes in the reservoir, indexes for Reservoir(i).MacroNodesID
    Reservoir(r).EntryNodesIndex = cell(1,NumModes); % border nodes for entry in the reservoir, indexes for Reservoir(i).MacroNodesID
    Reservoir(r).ExitNodesIndex = cell(1,NumModes); % border nodes for exit of the reservoir, indexes for Reservoir(i).MacroNodesID
    
    Reservoir(r).ModeIndex = cell(1,NumModes); % mode of each route passing through the reservoir 1) Car 2) Bus
    Reservoir(r).RouteMode = []; % mode of each route passing through the reservoir 1) Car 2) Bus
end

Temp_index = ones(1,NumRes);
for iroute = 1:NumRoutes
    i_p = 1;
    i_m = Route(iroute).ModeID;
    for r = Route(iroute).ResPath % loop on all reservoirs in the route
        Reservoir(r).TripLengthPerRoute = [Reservoir(r).TripLengthPerRoute Route(iroute).TripLengths(i_p)];
        Reservoir(r).FreeFlowTTPerRoute = [Reservoir(r).FreeFlowTTPerRoute Route(iroute).FreeFlowTTPerReservoir(i_p)];
        Reservoir(r).RoutesID = [Reservoir(r).RoutesID iroute];
        Reservoir(r).RoutesPathIndex = [Reservoir(r).RoutesPathIndex i_p];
        Reservoir(r).OriginRes = [Reservoir(r).OriginRes Route(iroute).ResOriginID];
        Reservoir(r).DestinationRes = [Reservoir(r).DestinationRes Route(iroute).ResDestinationID];
        Reservoir(r).ODmacroID{i_m} = [Reservoir(r).ODmacroID{i_m} Route(iroute).ODmacroID];
        
        inode1 = Route(iroute).NodePath(i_p); % entry node ID
        inode2 = Route(iroute).NodePath(i_p+1); % exit node ID
        i_n1 = find(Reservoir(r).MacroNodesID == inode1);
        i_n2 = find(Reservoir(r).MacroNodesID == inode2);
        Reservoir(r).RoutesNodeID = [Reservoir(r).RoutesNodeID [inode1; inode2]];
        Reservoir(r).NodeRoutesIndex{i_m,i_n1} = [Reservoir(r).NodeRoutesIndex{i_m,i_n1} Temp_index(r)];
        Reservoir(r).NodeRoutesIndex{i_m,i_n2} = [Reservoir(r).NodeRoutesIndex{i_m,i_n2} Temp_index(r)];
        
        if strcmp(MacroNode(inode1).Type,'origin')
            Reservoir(r).OriNodesIndex{i_m} = [Reservoir(r).OriNodesIndex{i_m} i_n1];
            Reservoir(r).OriRoutesIndex{i_m} = [Reservoir(r).OriRoutesIndex{i_m} Temp_index(r)];
        elseif strcmp(MacroNode(inode1).Type,'border') || strcmp(MacroNode(inode1).Type,'externalentry')
            Reservoir(r).EntryNodesIndex{i_m} = [Reservoir(r).EntryNodesIndex{i_m} i_n1];
            Reservoir(r).EntryRoutesIndex{i_m} = [Reservoir(r).EntryRoutesIndex{i_m} Temp_index(r)];
        end
        if strcmp(MacroNode(inode2).Type,'destination')
            Reservoir(r).DestNodesIndex{i_m} = [Reservoir(r).DestNodesIndex{i_m} i_n2];
            Reservoir(r).DestRoutesIndex{i_m} = [Reservoir(r).DestRoutesIndex{i_m} Temp_index(r)];
        elseif strcmp(MacroNode(inode2).Type,'border') || strcmp(MacroNode(inode2).Type,'externalexit')
            Reservoir(r).ExitNodesIndex{i_m} = [Reservoir(r).ExitNodesIndex{i_m} i_n2];
            Reservoir(r).ExitRoutesIndex{i_m} = [Reservoir(r).ExitRoutesIndex{i_m} Temp_index(r)];
        end
        if strcmp(Route(iroute).Mode, 'cartrip')
            Reservoir(r).ModeIndex{1} = [Reservoir(r).ModeIndex{1} Temp_index(r)];
        elseif strcmp(Route(iroute).Mode, 'publictransport')
            Reservoir(r).ModeIndex{2} = [Reservoir(r).ModeIndex{2} Temp_index(r)];
        end
        if strcmp(Route(iroute).Mode, 'cartrip')
            Reservoir(r).RouteMode = [Reservoir(r).RouteMode 1];
        elseif strcmp(Route(iroute).Mode, 'publictransport')
            Reservoir(r).RouteMode = [Reservoir(r).RouteMode 2];
        end
        Temp_index(r) = Temp_index(r) + 1;
        i_p = i_p + 1;
        
    end
end
for r = 1:NumRes
    for m = 1:NumModes
        if ~isempty(Reservoir(r).ODmacroID{m})
            Reservoir(r).ODmacroID{m} = unique(Reservoir(r).ODmacroID{m}); % remove repetitions
        end
        if ~isempty(Reservoir(r).OriNodesIndex{m})
            Reservoir(r).OriNodesIndex{m} = unique(Reservoir(r).OriNodesIndex{m});
        end
        if ~isempty(Reservoir(r).EntryNodesIndex{m})
            Reservoir(r).EntryNodesIndex{m} = unique(Reservoir(r).EntryNodesIndex{m});
        end
        if ~isempty(Reservoir(r).DestNodesIndex{m})
            Reservoir(r).DestNodesIndex{m} = unique(Reservoir(r).DestNodesIndex{m});
        end
        if ~isempty(Reservoir(r).ExitNodesIndex{m})
            Reservoir(r).ExitNodesIndex{m} = unique(Reservoir(r).ExitNodesIndex{m});
        end
    end
end
for iroute = 1:NumRoutes
    Route(iroute).ResRouteIndex = zeros(1,NumRes); % index of the route in the set of crossing routes in each reservoir
end
for r = 1:NumRes
    if ~isempty(Reservoir(r).RoutesID)
        i_r = 1;
        for iroute = Reservoir(r).RoutesID
            Route(iroute).ResRouteIndex(r) = i_r;
            i_r = i_r + 1;
        end
    end
end


% Variable initialization for convergence calculation
Assignment.OldAssignCoeff = zeros(1,NumRoutes);
Assignment.NewAssignCoeff = zeros(1,NumRoutes);
Assignment.OldMeanTravelTime = zeros(1,NumRoutes);
for iroute = 1:NumRoutes
    Route(iroute).ListAssignCoeff = [];
    Route(iroute).AssignCoeff = 0;
    Assignment.OldMeanTravelTime(iroute) = Route(iroute).TotalTime;
end


% Save structures for out files
Out.Reservoir = Reservoir;
Out.Route = Route;





