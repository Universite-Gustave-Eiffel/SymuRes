%% POST-PROCESSING
%--------------------------------------------------------------------------
% Accumulation-based model
%
% June 2020 - Guilhem Mariotte

NumRes = length(Reservoir);
NumODmacro = length(ODmacro);
NumRoutes = length(Route);
NumPaths = Assignment.NumShortestPath;
NumMacroNodes = length(MacroNode);

TimeStep = Simulation.TimeStep;
SimulTime = Simulation.Time;
NumTimes = floor(Simulation.Duration/TimeStep) + 1; % number of times
NumModes = Simulation.NumModes;

MFDfct = Simulation.MFDfct;


%% Deriving accumulation and flows for different levels of aggregation
%--------------------------------------------------------------------------

% Default null values
for r = 1:NumRes
    if isempty(Reservoir(r).RoutesID)
        Reservoir(r).AccPerRoute = zeros(1,NumTimes);
        Reservoir(r).InflowPerRoute = zeros(1,NumTimes);
        Reservoir(r).OutflowPerRoute = zeros(1,NumTimes);
        if isfield(Reservoir,'AccCircuPerRoute')
            Reservoir(r).AccCircuPerRoute = zeros(1,NumTimes);
            Reservoir(r).AccQueuePerRoute = zeros(1,NumTimes);
            Reservoir(r).OutflowCircuPerRoute = zeros(1,NumTimes);
        end
    end
end

% Variable initialization
for r = 1:NumRes
    Temp_Nroutes = max([length(Reservoir(r).RoutesID) 1]);
    
    Reservoir(r).Acc = zeros(1,NumTimes); % total accumulation [veh]
    Reservoir(r).AccPerMode = zeros(NumModes,NumTimes); % accumulation per mode [veh]
    Reservoir(r).MeanSpeed = zeros(1,NumTimes); % mean speed [m/s]
    Reservoir(r).MeanSpeedPerMode = zeros(NumModes,NumTimes); % mean speed per mode [m/s]
    Reservoir(r).AvgTripLength = zeros(1,NumTimes); % average trip length [m]
    Reservoir(r).InflowPerResPerDest = zeros(NumRes,NumRes,NumTimes); % (adj res, dest res, time)
    Reservoir(r).Inflow = zeros(1,NumTimes); % total effective inflow [veh/s]
    Reservoir(r).OutflowPerResPerDest = zeros(NumRes,NumRes,NumTimes); % (adj res, dest res, time)
    Reservoir(r).Outflow = zeros(1,NumTimes); % total effective outflow [veh/s]
    Reservoir(r).NinPerRoute = zeros(Temp_Nroutes,NumTimes); % entry cumulative count per route
    Reservoir(r).NinPerResPerDest = zeros(NumRes,NumRes,NumTimes); % nin(adj res,d,t)
    Reservoir(r).Nin = zeros(1,NumTimes); % entry cumulative count [veh]
    Reservoir(r).NoutPerRoute = zeros(Temp_Nroutes,NumTimes); % exit cumulative count per route
    Reservoir(r).NoutPerResPerDest = zeros(NumRes,NumRes,NumTimes); % nout(adj res,d,t)
    Reservoir(r).Nout = zeros(1,NumTimes); % exit cumulative count [veh]
    if isfield(Reservoir,'AccCircuPerRoute')
        Reservoir(r).NoutCircuPerRoute = zeros(Temp_Nroutes,NumTimes); % exit cumulative count of circulating veh per route
    end
end

% Total accumulation, inflow and outflow
for r = 1:NumRes
    % Total accumulation
    Reservoir(r).Acc = sum(Reservoir(r).AccPerRoute,1);
    
    % Accumulation per mode
    i_r = 1;
    for iroute = Reservoir(r).RoutesID
        i_m = Route(iroute).ModeID;
        Reservoir(r).AccPerMode(i_m,:) = Reservoir(r).AccPerMode(i_m,:) + Reservoir(r).AccPerRoute(i_r,:);
        i_r = i_r + 1;
    end
    
    % Mean speed
    Reservoir(r).MeanSpeedPerMode = Reservoir(r).MeanSpeed; % mean speed per mode [m/s]
    
    % Inflow/outflow
    Reservoir(r).Inflow = sum(Reservoir(r).InflowPerRoute,1);
    Reservoir(r).Outflow = sum(Reservoir(r).OutflowPerRoute,1);
end

% Inflow and outflow per adjacent reservoir
for r = 1:NumRes
    i_r = 1;
    for iroute = Reservoir(r).RoutesID % loop on all routes including Rr
        i = Route(iroute).ResOriginID;
        j = Route(iroute).ResDestinationID;
        
        Temp_path = Route(iroute).ResPath;
        if r == i && r == j % internal trip
            Temp_prevresID = i;
            Temp_nextresID = j;
        elseif r == i && r ~= j % if Rr is the route origin
            Temp_prevresID = 0;
            Temp_nextresID = Temp_path(Reservoir(r).RoutesPathIndex(i_r)+1); % next reservoir in the path
        elseif r == j && r ~= i % if Rr is the route destination
            Temp_prevresID = Temp_path(Reservoir(r).RoutesPathIndex(i_r)-1); % previous reservoir in the path
            Temp_nextresID = 0;
        else
            Temp_prevresID = Temp_path(Reservoir(r).RoutesPathIndex(i_r)-1); % previous reservoir in the path
            Temp_nextresID = Temp_path(Reservoir(r).RoutesPathIndex(i_r)+1); % next reservoir in the path
        end
        
        % Flows for adjacent reservoirs
        for iadj = Reservoir(r).AdjacentRes % loop on all adjacent reservoirs
            if iadj == Temp_prevresID && iadj == Temp_nextresID
                Temp_nin = Reservoir(r).AccPerRoute(i_r,1);
                Temp_qin = Reservoir(r).InflowPerRoute(i_r,:);
                Temp_qin2(1,1,:) = Temp_qin(1,:);
                Reservoir(r).InflowPerResPerDest(iadj,j,:) = Reservoir(r).InflowPerResPerDest(iadj,j,:) + Temp_qin2;
                Reservoir(r).NinPerResPerDest(iadj,j,1) = Temp_nin;
                Temp_qout = Reservoir(r).OutflowPerRoute(i_r,:);
                Temp_qout2(1,1,:) = Temp_qout(1,:);
                Reservoir(r).OutflowPerResPerDest(iadj,j,:) = Reservoir(r).OutflowPerResPerDest(iadj,j,:) + Temp_qout2;
            elseif iadj == Temp_prevresID % if the adjacent reservoir is included in the route
                Temp_nin = Reservoir(r).AccPerRoute(i_r,1);
                Temp_qin = Reservoir(r).InflowPerRoute(i_r,:);
                Temp_qin2(1,1,:) = Temp_qin(1,:);
                Reservoir(r).InflowPerResPerDest(iadj,j,:) = Reservoir(r).InflowPerResPerDest(iadj,j,:) + Temp_qin2;
                Reservoir(r).NinPerResPerDest(iadj,j,1) = Temp_nin;
            elseif iadj == Temp_nextresID
                Temp_qout = Reservoir(r).OutflowPerRoute(i_r,:);
                Temp_qout2(1,1,:) = Temp_qout(1,:);
                Reservoir(r).OutflowPerResPerDest(iadj,j,:) = Reservoir(r).OutflowPerResPerDest(iadj,j,:) + Temp_qout2;
            end
        end
        i_r = i_r + 1;
    end
end

% Cumulative count curves
for r = 1:NumRes
    Reservoir(r).NinPerRoute(:,1) = Reservoir(r).AccPerRoute(:,1);
    Reservoir(r).AvgTripLength(1) = mean(Reservoir(r).TripLengthPerRoute);
    
    for itime = 2:NumTimes
        Reservoir(r).NinPerRoute(:,itime) = Reservoir(r).NinPerRoute(:,itime-1) + TimeStep*Reservoir(r).InflowPerRoute(:,itime-1);
        Reservoir(r).NinPerResPerDest(:,:,itime) = Reservoir(r).NinPerResPerDest(:,:,itime-1) + TimeStep*Reservoir(r).InflowPerResPerDest(:,:,itime-1);
        Reservoir(r).Nin(itime) = Reservoir(r).Nin(itime-1) + TimeStep*Reservoir(r).Inflow(itime-1);
        
        Reservoir(r).NoutPerRoute(:,itime) = Reservoir(r).NoutPerRoute(:,itime-1) + TimeStep*Reservoir(r).OutflowPerRoute(:,itime-1);
        Reservoir(r).NoutPerResPerDest(:,:,itime) = Reservoir(r).NoutPerResPerDest(:,:,itime-1) + TimeStep*Reservoir(r).OutflowPerResPerDest(:,:,itime-1);
        Reservoir(r).Nout(itime) = Reservoir(r).Nout(itime-1) + TimeStep*Reservoir(r).Outflow(itime-1);
        
        if Reservoir(r).Acc(itime) > 0
            Reservoir(r).AvgTripLength(itime) = Reservoir(r).Acc(itime)/sum(Reservoir(r).AccPerRoute(:,itime)'./Reservoir(r).TripLengthPerRoute);
        else
            Reservoir(r).AvgTripLength(itime) = Reservoir(r).AvgTripLength(itime-1);
        end
    end
end

% ODmacro demand
for od = 1:NumODmacro
    for j = 1:length(ODmacro(od).Demand)
        Temp_times = ODmacro(od).Demand(j).Time;
        Temp_data = ODmacro(od).Demand(j).Data;
        [Temp_times, Temp_data] = stairfct(Temp_times,Temp_data,TimeStep,0,Simulation.Duration);
        ODmacro(od).Demand(j).Data2 = Temp_data; % demand sampled to the simulation timestep
    end
end

% Route demand
if Assignment.PredefRoute == 0
    for iroute = 1:NumRoutes
        od = Route(iroute).ODmacroID;
        for j = 1:length(ODmacro(od).Demand)
            if strcmp(ODmacro(od).Demand(j).Purpose,'cartrip')
                for iperiod = 1:(length(Assignment.Periods)-1)
                    Temp_StartTimeID = floor(Assignment.Periods(iperiod)/TimeStep) + 1;
                    Temp_EndTimeID = floor(Assignment.Periods(iperiod+1)/TimeStep) + 1;
                    Temp_demOD = ODmacro(od).Demand(j).Data2(Temp_StartTimeID:Temp_EndTimeID);
                    Route(iroute).Demand(Temp_StartTimeID:Temp_EndTimeID) = Route(iroute).ListAssignCoeff(iperiod).*Temp_demOD; % row vector
                end
            end
        end
    end
end

% for iroute = 1:NumRoutes
%     if Route(iroute).AssignCoeff > 0
%         o = Route(iroute).ResOriginID;
%         d = Route(iroute).ResDestinationID;
%         i_o = Route(iroute).ResRouteIndex(o);
%         i_d = Route(iroute).ResRouteIndex(d);
%         
%         Temp_nin = Reservoir(o).NinPerRoute(i_o,:);
%         Temp_nout = Reservoir(d).NoutPerRoute(i_d,:);
%         Route(iroute).TravelTime = ExperiencedTravelTime(CurrentTime,Temp_nin,Temp_nout);
%     else
%         % If the route is not used, default is the free-flow travel time
%         Route(iroute).TravelTime = Route(iroute).FreeFlowTravelTime*ones(1,NumTimes);
%     end
%     % TT for the last time step (simplification)
%     Route(iroute).TravelTime(NumTimes) = Route(iroute).TravelTime(NumTimes-1);
% end

% Dummy Vehicle structure
Vehicle = [];


clear Temp_*


