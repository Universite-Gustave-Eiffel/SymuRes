%% POST-PROCESSING
%--------------------------------------------------------------------------
% Trip-based model
%
% Nov 2019 - Guilhem Mariotte

NumRes = length(Reservoir);
NumODmacro = length(ODmacro);
NumRoutes = length(Route);
NumPaths = Assignment.NumShortestPath;
NumMacroNodes = length(MacroNode);
NumVeh = length(Vehicle);

TimeStep = Simulation.TimeStep;
SimulTime = Simulation.Time;
NumTimes = floor(Simulation.Duration/TimeStep) + 1; % number of times

MFDfct = Simulation.MFDfct;

eps0 = 0.5; % criterion for waiting times [s]


%% Deriving accumulation and N-curves from the Vehicle entry/exit times
%--------------------------------------------------------------------------

% Variable initialization
for r = 1:NumRes
    Temp_Nroutes = max([length(Reservoir(r).RoutesID) 1]);
    Reservoir(r).Acc = zeros(1,NumTimes); % total accumulation [veh]
    Reservoir(r).AccPerRoute = zeros(Temp_Nroutes,NumTimes); % accumulation per route crossing Rr [veh]
    Reservoir(r).MeanSpeed = zeros(1,NumTimes); % mean speed [m/s]
    Reservoir(r).InflowPerRoute = zeros(Temp_Nroutes,NumTimes); % inflow per route [veh/s]
    Reservoir(r).OutflowPerRoute = zeros(Temp_Nroutes,NumTimes); % outflow per route [veh/s]
    Reservoir(r).NinPerRoute = zeros(Temp_Nroutes,NumTimes); % entry cumulative count per route [veh]
    Reservoir(r).NoutPerRoute = zeros(Temp_Nroutes,NumTimes); % exit cumulative count per route [veh]
    
    Temp_ncircu = zeros(NumRes,NumTimes); % acc circu(r,t)
    Temp_sumavgtriplength = zeros(NumRes,NumTimes); % sum L(r,t)
    Temp_countavgtriplength = zeros(NumRes,NumTimes); % number of L(r,t)
end

% Restart the vehicle locations
for i = 1:NumVeh % loop on all vehicles
    ires = Route(Vehicle(i).RouteID).ResPath(1); % first reservoir in the route of veh i
    Vehicle(i).CurrentResID = ires;
    Vehicle(i).PathIndex = 1;
end

% Replay the simulation chronology
itime = 1;
for i = 1:Global.NumEvents % loop on all chronological events
    iveh = Global.VehID(i);
    iroute = Vehicle(iveh).RouteID;
    o = Route(iroute).ResOriginID;
    d = Route(iroute).ResDestinationID;
    ires = Vehicle(iveh).CurrentResID;
    i_r = Route(iroute).ResRouteIndex(ires);
    Temp_pathindex = Vehicle(iveh).PathIndex;
    
    if ires == o && o ~= d % the veh is in the first reservoir of the route: the current event is an entry
        while itime < NumTimes && SimulTime(itime) < Vehicle(iveh).EntryTimes(Temp_pathindex)
            for r = 1:NumRes
                Reservoir(r).AccPerRoute(:,itime+1) = Reservoir(r).AccPerRoute(:,itime);
                Reservoir(r).NinPerRoute(:,itime+1) = Reservoir(r).NinPerRoute(:,itime);
                Reservoir(r).NoutPerRoute(:,itime+1) = Reservoir(r).NoutPerRoute(:,itime);
                Temp_sumavgtriplength(r,itime+1) = Temp_sumavgtriplength(r,itime);
                Temp_countavgtriplength(r,itime+1) = Temp_countavgtriplength(r,itime);
            end
            itime = itime + 1;
        end
        Reservoir(ires).AccPerRoute(i_r,itime) = Reservoir(ires).AccPerRoute(i_r,itime) + 1;
        Reservoir(ires).NinPerRoute(i_r,itime) = Reservoir(ires).NinPerRoute(i_r,itime) + 1;
        Vehicle(iveh).CurrentResID = Route(iroute).ResPath(Temp_pathindex+1); % next reservoir in the path, indicate that the veh has already entered
        
%         if i < 2000
%             fprintf('%s%i \t %s%i \t %s%i \t %s%i \t %s%i \t %s%i \t %s%i \t %s%3.1f \t %s%3.1f \n',...
%                 'case',1,'i=',i,'iveh=',iveh,'r=',ires,'o=',o,'d=',d,'p=',p,'t=',CurrentTime(itime),'n=',Reservoir(ires).AccPerODPerPath(o,d,p,itime))
%         end
        
    elseif ires == d && o ~= d && Temp_pathindex > 1 % the veh is in the last reservoir of the route: the current event is an exit
        while itime < NumTimes && SimulTime(itime) < Vehicle(iveh).ExitTimes(Temp_pathindex)
            for r = 1:NumRes
                Reservoir(r).AccPerRoute(:,itime+1) = Reservoir(r).AccPerRoute(:,itime);
                Reservoir(r).NinPerRoute(:,itime+1) = Reservoir(r).NinPerRoute(:,itime);
                Reservoir(r).NoutPerRoute(:,itime+1) = Reservoir(r).NoutPerRoute(:,itime);
                Temp_sumavgtriplength(r,itime+1) = Temp_sumavgtriplength(r,itime);
                Temp_countavgtriplength(r,itime+1) = Temp_countavgtriplength(r,itime);
            end
            itime = itime + 1;
        end
        Reservoir(ires).AccPerRoute(i_r,itime) = Reservoir(ires).AccPerRoute(i_r,itime) - 1;
        Reservoir(ires).NoutPerRoute(i_r,itime) = Reservoir(ires).NoutPerRoute(i_r,itime) + 1;
        Temp_sumavgtriplength(ires,itime) = Temp_sumavgtriplength(ires,itime) + Vehicle(iveh).TripLength(Temp_pathindex);
        Temp_countavgtriplength(ires,itime) = Temp_countavgtriplength(ires,itime) + 1;
        
    elseif o == d % if the route is actually one reservoir (internal trip): the current event may be an entry or an exit
        if Temp_pathindex == 1 % entry in the reservoir (creation)
            while itime < NumTimes && SimulTime(itime) < Vehicle(iveh).EntryTimes(Temp_pathindex)
                for r = 1:NumRes
                    Reservoir(r).AccPerRoute(:,itime+1) = Reservoir(r).AccPerRoute(:,itime);
                    Reservoir(r).NinPerRoute(:,itime+1) = Reservoir(r).NinPerRoute(:,itime);
                    Reservoir(r).NoutPerRoute(:,itime+1) = Reservoir(r).NoutPerRoute(:,itime);
                    Temp_sumavgtriplength(r,itime+1) = Temp_sumavgtriplength(r,itime);
                    Temp_countavgtriplength(r,itime+1) = Temp_countavgtriplength(r,itime);
                end
                itime = itime + 1;
            end
            Reservoir(ires).AccPerRoute(i_r,itime) = Reservoir(ires).AccPerRoute(i_r,itime) + 1;
            Reservoir(ires).NinPerRoute(i_r,itime) = Reservoir(ires).NinPerRoute(i_r,itime) + 1;
            Vehicle(iveh).PathIndex = Temp_pathindex + 1; % indicate that the veh has already entered
            
        elseif Temp_pathindex == 2 % indicate the exit from the reservoir
            Temp_pathindex = 1; % correct value of the index
            while itime < NumTimes && SimulTime(itime) < Vehicle(iveh).ExitTimes(Temp_pathindex)
                for r = 1:NumRes
                    Reservoir(r).AccPerRoute(:,itime+1) = Reservoir(r).AccPerRoute(:,itime);
                    Reservoir(r).NinPerRoute(:,itime+1) = Reservoir(r).NinPerRoute(:,itime);
                    Reservoir(r).NoutPerRoute(:,itime+1) = Reservoir(r).NoutPerRoute(:,itime);
                    Temp_sumavgtriplength(r,itime+1) = Temp_sumavgtriplength(r,itime);
                    Temp_countavgtriplength(r,itime+1) = Temp_countavgtriplength(r,itime);
                end
                itime = itime + 1;
            end
            Reservoir(ires).AccPerRoute(i_r,itime) = Reservoir(ires).AccPerRoute(i_r,itime) - 1;
            Reservoir(ires).NoutPerRoute(i_r,itime) = Reservoir(ires).NoutPerRoute(i_r,itime) + 1;
            Temp_sumavgtriplength(ires,itime) = Temp_sumavgtriplength(ires,itime) + Vehicle(iveh).TripLength(Temp_pathindex);
            Temp_countavgtriplength(ires,itime) = Temp_countavgtriplength(ires,itime) + 1;
        end
        
    else % if the veh is in the middle of a multi-reservoir route: the current event is an exit and an entry
        if Temp_pathindex == 1 % if the veh is actually in the origin reservoir
            ires = o; % correct ID of the current reservoir
            i_r = Route(iroute).ResRouteIndex(ires);
        end
        while itime < NumTimes && SimulTime(itime) < Vehicle(iveh).ExitTimes(Temp_pathindex)
            for r = 1:NumRes
                Reservoir(r).AccPerRoute(:,itime+1) = Reservoir(r).AccPerRoute(:,itime);
                Reservoir(r).NinPerRoute(:,itime+1) = Reservoir(r).NinPerRoute(:,itime);
                Reservoir(r).NoutPerRoute(:,itime+1) = Reservoir(r).NoutPerRoute(:,itime);
                Temp_sumavgtriplength(r,itime+1) = Temp_sumavgtriplength(r,itime);
                Temp_countavgtriplength(r,itime+1) = Temp_countavgtriplength(r,itime);
            end
            itime = itime + 1;
        end
        Reservoir(ires).AccPerRoute(i_r,itime) = Reservoir(ires).AccPerRoute(i_r,itime) - 1;
        Reservoir(ires).NoutPerRoute(i_r,itime) = Reservoir(ires).NoutPerRoute(i_r,itime) + 1;
        Temp_sumavgtriplength(ires,itime) = Temp_sumavgtriplength(ires,itime) + Vehicle(iveh).TripLength(Temp_pathindex);
        Temp_countavgtriplength(ires,itime) = Temp_countavgtriplength(ires,itime) + 1;
        
        Vehicle(iveh).CurrentResID = Route(iroute).ResPath(Temp_pathindex+1); % next reservoir in the path
        Vehicle(iveh).PathIndex = Temp_pathindex + 1; % next reservoir in the path
        ires = Vehicle(iveh).CurrentResID;
        i_r = Route(iroute).ResRouteIndex(ires);
        Reservoir(ires).AccPerRoute(i_r,itime) = Reservoir(ires).AccPerRoute(i_r,itime) + 1;
        Reservoir(ires).NinPerRoute(i_r,itime) = Reservoir(ires).NinPerRoute(i_r,itime) + 1;
    end
%     if i < 2000
%         fprintf('%s%3.1f \t %s%i \t %s%i \t %s%i \t %s%i \t %s%i \t %s%i \t %s%i \t %s%i \n',...
%             't=',Global.SimulTime(i),'itime=',itime,'r=',ires,'o=',o,'d=',d,'p=',p,'n=',Reservoir(ires).AccPerODPerPath(o,d,p,itime),'sumL=',Temp_sumavgtriplength(ires,itime),'nbL=',Temp_countavgtriplength(ires,itime))
%     end
end

if itime < NumTimes
    for i = itime:(NumTimes-1)
        for r = 1:NumRes
            Reservoir(r).AccPerRoute(:,i+1) = Reservoir(r).AccPerRoute(:,i);
            Reservoir(r).NinPerRoute(:,i+1) = Reservoir(r).NinPerRoute(:,i);
            Reservoir(r).NoutPerRoute(:,i+1) = Reservoir(r).NoutPerRoute(:,i);
            Temp_sumavgtriplength(r,i+1) = Temp_sumavgtriplength(r,i);
            Temp_countavgtriplength(r,i+1) = Temp_countavgtriplength(r,i);
        end
    end
end

% Apply the simulation scaling factor
for r = 1:NumRes
    Temp_Nroutes = length(Reservoir(r).RoutesID);
    Reservoir(r).Acc = 1./Simulation.TripbasedSimuFactor.*Reservoir(r).Acc;
    Reservoir(r).AccPerRoute = 1./Simulation.TripbasedSimuFactor.*Reservoir(r).AccPerRoute;
    Reservoir(r).InflowPerRoute = 1./Simulation.TripbasedSimuFactor.*Reservoir(r).InflowPerRoute;
    Reservoir(r).OutflowPerRoute = 1./Simulation.TripbasedSimuFactor.*Reservoir(r).OutflowPerRoute;
    Reservoir(r).NinPerRoute = 1./Simulation.TripbasedSimuFactor.*Reservoir(r).NinPerRoute;
    Reservoir(r).NoutPerRoute = 1./Simulation.TripbasedSimuFactor.*Reservoir(r).NoutPerRoute;
    
    Temp_countavgtriplength = 1./Simulation.TripbasedSimuFactor.*Temp_countavgtriplength;
end



%% Deriving accumulation and flows for different levels of aggregation
%--------------------------------------------------------------------------

% Variable initialization
for r = 1:NumRes
    Reservoir(r).InflowPerResPerDest = zeros(NumRes,NumRes,NumTimes); % (adj res, dest res, time)
    Reservoir(r).Inflow = zeros(1,NumTimes); % total effective inflow
    Reservoir(r).OutflowPerResPerDest = zeros(NumRes,NumRes,NumTimes); % (adj res, dest res, time)
    Reservoir(r).Outflow = zeros(1,NumTimes); % total effective outflow
    Reservoir(r).NinPerResPerDest = zeros(NumRes,NumRes,NumTimes); % nin(adj res,d,t)
    Reservoir(r).Nin = zeros(1,NumTimes); % nin(t)
    Reservoir(r).NoutPerResPerDest = zeros(NumRes,NumRes,NumTimes); % nout(adj res,d,t)
    Reservoir(r).Nout = zeros(1,NumTimes); % nout(t)
end

% Inflow and outflow per route
Ni = 10; % adopt a different sampling to avoid too many oscillations of flows
% Ni = 5;
ti = SimulTime(1:Ni:NumTimes);
for r = 1:NumRes
    for i_r = 1:length(Reservoir(r).RoutesID)
        Temp_inflow = deriv(ti,Reservoir(r).NinPerRoute(i_r,1:Ni:NumTimes));
        Temp_outflow = deriv(ti,Reservoir(r).NoutPerRoute(i_r,1:Ni:NumTimes));
        Reservoir(r).InflowPerRoute(i_r,:) = resamp(SimulTime,ti,Temp_inflow);
        Reservoir(r).OutflowPerRoute(i_r,:) = resamp(SimulTime,ti,Temp_outflow);
        %Reservoir(r).InflowPerRoute(i_r,:) = ExpoMovingAvg(Reservoir(r).InflowPerRoute(i_r,:),0.01);
        %Reservoir(r).OutflowPerRoute(i_r,:) = ExpoMovingAvg(Reservoir(r).OutflowPerRoute(i_r,:),0.01);
    end
end

% Circulating accumulation per reservoir
% for itime = 1:NumTimes
%     for iveh = 1:NumVeh
%         Temp_Npath = length(Vehicle(iveh).TripLength);
%         for i_p = 1:Temp_Npath % loop on all res in the path
%             % if the veh has not exited yet
%             if Vehicle(iveh).EntryTimes(i_p) <= CurrentTime(itime) && ...
%                     Vehicle(iveh).ExitTimes(i_p) == Inf
%                 if Vehicle(iveh).WaitingTimes(i_p) > 0 % if the veh is waiting for exiting
%                     ires = Route(Vehicle(iveh).RouteID).Path(i_p);
%                     Temp_ncircu(ires,itime) = Temp_ncircu(ires,itime) - 1;
%                 end
%             else
%                 if Vehicle(iveh).ExitTimes(i_p)-Vehicle(iveh).WaitingTimes(i_p) <= CurrentTime(itime) && ...
%                         CurrentTime(itime) < Vehicle(iveh).ExitTimes(i_p) % if the veh is waiting for exiting
%                     ires = Route(Vehicle(iveh).RouteID).Path(i_p);
%                     Temp_ncircu(ires,itime) = Temp_ncircu(ires,itime) - 1;
%                 end
%             end
%         end
%     end
% end

% Temp_vehlist = [];
% iveh0 = 0;
% for itime = 1:NumTimes
%     % Append newly created vehicles in the list
%     iveh0 = iveh0 + 1;
%     while iveh0 <= NumVeh && Vehicle(iveh0).CreationTime <= CurrentTime(itime)
%         Temp_vehlist = [Temp_vehlist iveh0];
%         iveh0 = iveh0 + 1;
%     end
%     iveh0 = iveh0 - 1;
%     % Remove exited vehicles from the list
%     if iveh0 > 0
%         Temp_indlist = [];
%         i = 1;
%         for iveh = Temp_vehlist
%             if Vehicle(iveh).ExitTimes(end) >= CurrentTime(itime)
%                 Temp_indlist = [Temp_indlist i];
%             end
%             i = i + 1;
%         end
%         Temp_vehlist = Temp_vehlist(Temp_indlist);
%     end
%     %disp(int2str(length(Temp_vehlist)))
%     % Look for waiting vehicles
%     for iveh = Temp_vehlist
%         Temp_Npath = length(Vehicle(iveh).TripLength);
%         for i_p = 1:Temp_Npath % loop on all res in the path
%             % if the veh has not exited yet
%             if Vehicle(iveh).EntryTimes(i_p) <= CurrentTime(itime) && ...
%                     Vehicle(iveh).ExitTimes(i_p) == Inf
%                 if Vehicle(iveh).WaitingTimes(i_p) > 0 % if the veh is waiting for exiting
%                     ires = Route(Vehicle(iveh).RouteID).Path(i_p);
%                     Temp_ncircu(ires,itime) = Temp_ncircu(ires,itime) - 1;
%                 end
%             else
%                 if Vehicle(iveh).ExitTimes(i_p)-Vehicle(iveh).WaitingTimes(i_p) <= CurrentTime(itime) && ...
%                         CurrentTime(itime) < Vehicle(iveh).ExitTimes(i_p) % if the veh is waiting for exiting
%                     ires = Route(Vehicle(iveh).RouteID).Path(i_p);
%                     Temp_ncircu(ires,itime) = Temp_ncircu(ires,itime) - 1;
%                 end
%             end
%         end
%     end
% end

% Total accumulation and mean speed
for r = 1:NumRes
    Reservoir(r).Acc = sum(Reservoir(r).AccPerRoute,1);
    Temp_param = Reservoir(r).MFDfctParam;
    for i = 1:NumTimes
        Temp_nr = Reservoir(r).Acc(i);
        %Temp_nr = Temp_ncircu(r,i);
        if Temp_nr == 0
            Temp_Vr = Reservoir(r).FreeflowSpeed;
        else
            Temp_Vr = MFDfct(Temp_nr,Temp_param)/Temp_nr;
        end
        Reservoir(r).MeanSpeed(i) = Temp_Vr;
    end
end

% Average trip lengths
for r = 1:NumRes
    Reservoir(r).AvgTripLength = Temp_sumavgtriplength(r,:)./Temp_countavgtriplength(r,:);
end

% Total inflow and outflow
for r = 1:NumRes
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
    for itime = 2:NumTimes
        Reservoir(r).NinPerResPerDest(:,:,itime) = Reservoir(r).NinPerResPerDest(:,:,itime-1) + TimeStep*Reservoir(r).InflowPerResPerDest(:,:,itime-1);
        Reservoir(r).Nin(itime) = Reservoir(r).Nin(itime-1) + TimeStep*Reservoir(r).Inflow(itime-1);
        
        Reservoir(r).NoutPerResPerDest(:,:,itime) = Reservoir(r).NoutPerResPerDest(:,:,itime-1) + TimeStep*Reservoir(r).OutflowPerResPerDest(:,:,itime-1);
        Reservoir(r).Nout(itime) = Reservoir(r).Nout(itime-1) + TimeStep*Reservoir(r).Outflow(itime-1);
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


clear Temp_*

