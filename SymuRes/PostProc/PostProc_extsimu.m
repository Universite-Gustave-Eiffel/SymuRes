%% POST-PROCESSING
%--------------------------------------------------------------------------

NumRes = length(Reservoir);
NumODmacro = length(ODmacro);
NumRoutes = length(Route);
NumPaths = 3;
NumVeh = length(Vehicle);

SimulationDuration = Simulation.Duration;
TimeStep = Simulation.TimeStep;
SimulTime = Simulation.Time;
NumTimes = floor(SimulationDuration/TimeStep) + 1; % number of times

MFDfct = Simulation.MFDfct;

eps0 = 0.5; % criterion for waiting times [s]

%% Deriving accumulation and N-curves from the Vehicle entry/exit times
%--------------------------------------------------------------------------

% Variable initialization
for r = 1:NumRes
    Reservoir(r).AccPerODPerPath = zeros(NumRes,NumRes,NumPaths,NumTimes); % accumulation per OD per path: n(o,d,p,t)
    Reservoir(r).NinPerODPerPath = zeros(NumRes,NumRes,NumPaths,NumTimes); % nin(o,d,p,t)
    Reservoir(r).NoutPerODPerPath = zeros(NumRes,NumRes,NumPaths,NumTimes); % nout(o,d,p,t)
    
    Reservoir(r).AccPerODPerPath(:,:,:,1) = 0; % initial accumulation
    Reservoir(r).NinPerODPerPath(:,:,:,1) = 0;
    
    Temp_ncircu = zeros(NumRes,NumTimes); % acc circu(r,t)
    Temp_sumtriplengthperOD = zeros(NumRes,NumRes,NumRes,NumTimes); % sum L(r,o,d,t)
    Temp_sumtriplengthperdest = zeros(NumRes,NumRes,NumTimes); % sum L(r,d,t)
    Temp_sumavgtriplength = zeros(NumRes,NumTimes); % sum L(r,t)
    Temp_counttriplengthperOD = zeros(NumRes,NumRes,NumRes,NumTimes); % number of L(r,o,d,t)
    Temp_counttriplengthperdest = zeros(NumRes,NumRes,NumTimes); % number of L(r,d,t)
    Temp_countavgtriplength = zeros(NumRes,NumTimes); % number of L(r,t)
end

% Restart the vehicle locations
for i = 1:NumVeh % loop on all vehicles
    ires = Route(Vehicle(i).RouteID).Path(1); % first reservoir in the route of veh i
    Vehicle(i).CurrentResID = ires;
    Vehicle(i).PathIndex = 1;
end

% Replay the simulation chronology
itime = 1;
for i = 1:Global.NumEvents % loop on all chronological events
    iveh = Global.VehID(i);
    iroute = Vehicle(iveh).RouteID;
    o = Route(iroute).OriginID;
    d = Route(iroute).DestinationID;
    p = Route(iroute).PathID;
    ires = Vehicle(iveh).CurrentResID;
    Temp_pathindex = Vehicle(iveh).PathIndex;
    
    if ires == o && o ~= d % the veh is in the first reservoir of the route: the current event is an entry
        while itime < NumTimes && SimulTime(itime) < Vehicle(iveh).EntryTimes(Temp_pathindex)
            for r = 1:NumRes
                Reservoir(r).AccPerODPerPath(:,:,:,itime+1) = Reservoir(r).AccPerODPerPath(:,:,:,itime);
                Reservoir(r).NinPerODPerPath(:,:,:,itime+1) = Reservoir(r).NinPerODPerPath(:,:,:,itime);
                Reservoir(r).NoutPerODPerPath(:,:,:,itime+1) = Reservoir(r).NoutPerODPerPath(:,:,:,itime);
                Temp_sumavgtriplength(r,itime+1) = Temp_sumavgtriplength(r,itime);
                Temp_countavgtriplength(r,itime+1) = Temp_countavgtriplength(r,itime);
            end
            itime = itime + 1;
        end
        Reservoir(ires).AccPerODPerPath(o,d,p,itime) = Reservoir(ires).AccPerODPerPath(o,d,p,itime) + 1;
        Reservoir(ires).NinPerODPerPath(o,d,p,itime) = Reservoir(ires).NinPerODPerPath(o,d,p,itime) + 1;
        Vehicle(iveh).CurrentResID = Route(iroute).Path(Temp_pathindex+1); % next reservoir in the path, indicate that the veh has already entered
        
%         if i < 2000
%             fprintf('%s%i \t %s%i \t %s%i \t %s%i \t %s%i \t %s%i \t %s%i \t %s%3.1f \t %s%3.1f \n',...
%                 'case',1,'i=',i,'iveh=',iveh,'r=',ires,'o=',o,'d=',d,'p=',p,'t=',CurrentTime(itime),'n=',Reservoir(ires).AccPerODPerPath(o,d,p,itime))
%         end
        
    elseif ires == d && o ~= d && Temp_pathindex > 1 % the veh is in the last reservoir of the route: the current event is an exit
        while itime < NumTimes && SimulTime(itime) < Vehicle(iveh).ExitTimes(Temp_pathindex)
            for r = 1:NumRes
                Reservoir(r).AccPerODPerPath(:,:,:,itime+1) = Reservoir(r).AccPerODPerPath(:,:,:,itime);
                Reservoir(r).NinPerODPerPath(:,:,:,itime+1) = Reservoir(r).NinPerODPerPath(:,:,:,itime);
                Reservoir(r).NoutPerODPerPath(:,:,:,itime+1) = Reservoir(r).NoutPerODPerPath(:,:,:,itime);
                Temp_sumavgtriplength(r,itime+1) = Temp_sumavgtriplength(r,itime);
                Temp_countavgtriplength(r,itime+1) = Temp_countavgtriplength(r,itime);
            end
            itime = itime + 1;
        end
        Reservoir(ires).AccPerODPerPath(o,d,p,itime) = Reservoir(ires).AccPerODPerPath(o,d,p,itime) - 1;
        Reservoir(ires).NoutPerODPerPath(o,d,p,itime) = Reservoir(ires).NoutPerODPerPath(o,d,p,itime) + 1;
        Temp_sumavgtriplength(ires,itime) = Temp_sumavgtriplength(ires,itime) + Vehicle(iveh).TripLength(Temp_pathindex);
        Temp_countavgtriplength(ires,itime) = Temp_countavgtriplength(ires,itime) + 1;
        
    elseif o == d % if the route is actually one reservoir (internal trip): the current event may be an entry or an exit
        if Temp_pathindex == 1 % entry in the reservoir (creation)
            while itime < NumTimes && SimulTime(itime) < Vehicle(iveh).EntryTimes(Temp_pathindex)
                for r = 1:NumRes
                    Reservoir(r).AccPerODPerPath(:,:,:,itime+1) = Reservoir(r).AccPerODPerPath(:,:,:,itime);
                    Reservoir(r).NinPerODPerPath(:,:,:,itime+1) = Reservoir(r).NinPerODPerPath(:,:,:,itime);
                    Reservoir(r).NoutPerODPerPath(:,:,:,itime+1) = Reservoir(r).NoutPerODPerPath(:,:,:,itime);
                    Temp_sumavgtriplength(r,itime+1) = Temp_sumavgtriplength(r,itime);
                    Temp_countavgtriplength(r,itime+1) = Temp_countavgtriplength(r,itime);
                end
                itime = itime + 1;
            end
            Reservoir(ires).AccPerODPerPath(o,d,p,itime) = Reservoir(ires).AccPerODPerPath(o,d,p,itime) + 1;
            Reservoir(ires).NinPerODPerPath(o,d,p,itime) = Reservoir(ires).NinPerODPerPath(o,d,p,itime) + 1;
            Vehicle(iveh).PathIndex = Temp_pathindex + 1; % indicate that the veh has already entered
            
        elseif Temp_pathindex == 2 % indicate the exit from the reservoir
            Temp_pathindex = 1; % correct value of the index
            while itime < NumTimes && SimulTime(itime) < Vehicle(iveh).ExitTimes(Temp_pathindex)
                for r = 1:NumRes
                    Reservoir(r).AccPerODPerPath(:,:,:,itime+1) = Reservoir(r).AccPerODPerPath(:,:,:,itime);
                    Reservoir(r).NinPerODPerPath(:,:,:,itime+1) = Reservoir(r).NinPerODPerPath(:,:,:,itime);
                    Reservoir(r).NoutPerODPerPath(:,:,:,itime+1) = Reservoir(r).NoutPerODPerPath(:,:,:,itime);
                    Temp_sumavgtriplength(r,itime+1) = Temp_sumavgtriplength(r,itime);
                    Temp_countavgtriplength(r,itime+1) = Temp_countavgtriplength(r,itime);
                end
                itime = itime + 1;
            end
            Reservoir(ires).AccPerODPerPath(o,d,p,itime) = Reservoir(ires).AccPerODPerPath(o,d,p,itime) - 1;
            Reservoir(ires).NoutPerODPerPath(o,d,p,itime) = Reservoir(ires).NoutPerODPerPath(o,d,p,itime) + 1;
            Temp_sumavgtriplength(ires,itime) = Temp_sumavgtriplength(ires,itime) + Vehicle(iveh).TripLength(Temp_pathindex);
            Temp_countavgtriplength(ires,itime) = Temp_countavgtriplength(ires,itime) + 1;
        end
        
    else % if the veh is in the middle of a multi-reservoir route: the current event is an exit and an entry
        if Temp_pathindex == 1 % if the veh is actually in the origin reservoir
            ires = o; % correct ID of the current reservoir
        end
        while itime < NumTimes && SimulTime(itime) < Vehicle(iveh).ExitTimes(Temp_pathindex)
            for r = 1:NumRes
                Reservoir(r).AccPerODPerPath(:,:,:,itime+1) = Reservoir(r).AccPerODPerPath(:,:,:,itime);
                Reservoir(r).NinPerODPerPath(:,:,:,itime+1) = Reservoir(r).NinPerODPerPath(:,:,:,itime);
                Reservoir(r).NoutPerODPerPath(:,:,:,itime+1) = Reservoir(r).NoutPerODPerPath(:,:,:,itime);
                Temp_sumavgtriplength(r,itime+1) = Temp_sumavgtriplength(r,itime);
                Temp_countavgtriplength(r,itime+1) = Temp_countavgtriplength(r,itime);
            end
            itime = itime + 1;
        end
        Reservoir(ires).AccPerODPerPath(o,d,p,itime) = Reservoir(ires).AccPerODPerPath(o,d,p,itime) - 1;
        Reservoir(ires).NoutPerODPerPath(o,d,p,itime) = Reservoir(ires).NoutPerODPerPath(o,d,p,itime) + 1;
        Temp_sumavgtriplength(ires,itime) = Temp_sumavgtriplength(ires,itime) + Vehicle(iveh).TripLength(Temp_pathindex);
        Temp_countavgtriplength(ires,itime) = Temp_countavgtriplength(ires,itime) + 1;
        
        Vehicle(iveh).CurrentResID = Route(iroute).Path(Temp_pathindex+1); % next reservoir in the path
        Vehicle(iveh).PathIndex = Temp_pathindex + 1; % next reservoir in the path
        ires = Vehicle(iveh).CurrentResID;
        Reservoir(ires).AccPerODPerPath(o,d,p,itime) = Reservoir(ires).AccPerODPerPath(o,d,p,itime) + 1;
        Reservoir(ires).NinPerODPerPath(o,d,p,itime) = Reservoir(ires).NinPerODPerPath(o,d,p,itime) + 1;
    end
%     if i < 2000
%         fprintf('%s%3.1f \t %s%i \t %s%i \t %s%i \t %s%i \t %s%i \t %s%i \t %s%i \t %s%i \n',...
%             't=',Global.SimulTime(i),'itime=',itime,'r=',ires,'o=',o,'d=',d,'p=',p,'n=',Reservoir(ires).AccPerODPerPath(o,d,p,itime),'sumL=',Temp_sumavgtriplength(ires,itime),'nbL=',Temp_countavgtriplength(ires,itime))
%     end
end

if itime < NumTimes
    for i = itime:(NumTimes-1)
        for r = 1:NumRes
            Reservoir(r).AccPerODPerPath(:,:,:,i+1) = Reservoir(r).AccPerODPerPath(:,:,:,i);
            Reservoir(r).NinPerODPerPath(:,:,:,i+1) = Reservoir(r).NinPerODPerPath(:,:,:,i);
            Reservoir(r).NoutPerODPerPath(:,:,:,i+1) = Reservoir(r).NoutPerODPerPath(:,:,:,i);
            Temp_sumavgtriplength(r,i+1) = Temp_sumavgtriplength(r,i);
            Temp_countavgtriplength(r,i+1) = Temp_countavgtriplength(r,i);
        end
    end
end



%% Deriving accumulation and flows for different levels of aggregation
%--------------------------------------------------------------------------

% Variable initialization
for r = 1:NumRes
    Reservoir(r).Acc = zeros(1,NumTimes); % total accumulation: n(t)
    Reservoir(r).AccPerDest = zeros(NumRes,NumTimes); % accumulation per destination: n(d,t)
    Reservoir(r).AccPerOD = zeros(NumRes,NumRes,NumTimes); % accumulation per destination: n(o,d,t)
    Reservoir(r).MeanSpeed = zeros(1,NumTimes); % mean speed [m/s]: V(t)
    
    Reservoir(r).InflowPerODPerPath = zeros(NumRes,NumRes,NumPaths,NumTimes); % qin(o,d,p,t)
    Reservoir(r).InflowPerResPerDest = zeros(NumRes,NumRes,NumTimes); % (adj res, dest res, time)
    Reservoir(r).InflowInternalPerDest = zeros(NumRes,NumTimes); % effective internal inflow per destination
    Reservoir(r).InflowTransferPerDest = zeros(NumRes,NumTimes);
    Reservoir(r).Inflow = zeros(1,NumTimes); % total effective inflow
    
    Reservoir(r).OutflowPerODPerPath = zeros(NumRes,NumRes,NumPaths,NumTimes); % qout(o,d,p,t)
    Reservoir(r).OutflowPerResPerDest = zeros(NumRes,NumRes,NumTimes); % (adj res, dest res, time)
    Reservoir(r).OutflowPerDest = zeros(NumRes,NumTimes);
    Reservoir(r).Outflow = zeros(1,NumTimes); % total effective outflow
    
    Reservoir(r).NinPerResPerDest = zeros(NumRes,NumRes,NumTimes); % nin(adj res,d,t)
    Reservoir(r).NinInternalPerDest = zeros(NumRes,NumTimes); % nin(d,t)
    Reservoir(r).NinTransferPerDest = zeros(NumRes,NumTimes); % nin(d,t)
    Reservoir(r).Nin = zeros(1,NumTimes); % nin(t)
    
    Reservoir(r).NoutPerResPerDest = zeros(NumRes,NumRes,NumTimes); % nout(adj res,d,t)
    Reservoir(r).NoutPerDest = zeros(NumRes,NumTimes); % nout(d,t)
    Reservoir(r).Nout = zeros(1,NumTimes); % nout(t)
end

% Inflow and outflow per OD per path
Ni = 10; % adopt a different sampling to avoid too many oscillations of flows
% Ni = 5;
ti = SimulTime(1:Ni:NumTimes);
for r = 1:NumRes
    for j = 1:NumRes
        for i = 1:NumRes
            for k = 1:NumPaths
                Temp_inflow = deriv(ti,Reservoir(r).NinPerODPerPath(i,j,k,1:Ni:NumTimes));
                Temp_outflow = deriv(ti,Reservoir(r).NoutPerODPerPath(i,j,k,1:Ni:NumTimes));
                Reservoir(r).InflowPerODPerPath(i,j,k,:) = resamp(SimulTime,ti,Temp_inflow);
                Reservoir(r).OutflowPerODPerPath(i,j,k,:) = resamp(SimulTime,ti,Temp_outflow);
                %Reservoir(r).InflowPerODPerPath(i,j,k,:) = ExpoMovingAvg(Reservoir(r).InflowPerODPerPath(i,j,k,:),0.01);
                %Reservoir(r).OutflowPerODPerPath(i,j,k,:) = ExpoMovingAvg(Reservoir(r).OutflowPerODPerPath(i,j,k,:),0.01);
            end
        end
    end
end

% Accumulation and inflow/outflow
for r = 1:NumRes
    Reservoir(r).AccPerOD(:,:,:) = sum(Reservoir(r).AccPerODPerPath,3);
	Reservoir(r).AccPerDest(:,:) = sum(Reservoir(r).AccPerOD,1);
    Reservoir(r).Acc = sum(Reservoir(r).AccPerDest,1);
    Temp_ncircu(r,:) = Reservoir(r).Acc;
	
	Reservoir(r).OutflowPerDest(:,:) = sum(sum(Reservoir(r).OutflowPerODPerPath,3),1);
    Reservoir(r).Outflow(1,:) = sum(Reservoir(r).OutflowPerDest,1);
    
	Reservoir(r).InflowInternalPerDest(:,:) = sum(Reservoir(r).InflowPerODPerPath(r,:,:,:),3);
	
	[i_, Temp_res] = find(1:NumRes ~= r);
    for i = Temp_res % loop on all reservoirs except from Rr
	    Temp_qin = zeros(NumRes,NumTimes);
	    Temp_qin(:,:) = sum(Reservoir(r).InflowPerODPerPath(i,:,:,:),3);
        Reservoir(r).InflowTransferPerDest(:,:) = Reservoir(r).InflowTransferPerDest(:,:) + Temp_qin;
    end
	
	Reservoir(r).Inflow(1,:) = sum(Reservoir(r).InflowInternalPerDest,1) + sum(Reservoir(r).InflowTransferPerDest,1);
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

% Mean speed
% for r = 1:NumRes
%     Temp_nj = Reservoir(r).MaxAcc;
%     Temp_nc = Reservoir(r).CritAcc;
%     Temp_Pc = Reservoir(r).MaxProd;
%     for i = 1:NumTimes
%         Temp_nr = Reservoir(r).Acc(i);
%         %Temp_nr = Temp_ncircu(r,i);
%         if Temp_nr == 0
%             Temp_Vr = Reservoir(r).FreeflowSpeed;
%         else
%             Temp_Vr = MFDfct(Temp_nr,Temp_nj,Temp_nc,Temp_Pc)/Temp_nr;
%         end
%         Reservoir(r).MeanSpeed(i) = Temp_Vr;
%     end
% end

% Inflow and outflow per adjacent reservoir
for r = 1:NumRes
    Temp_index = 1;
    for iroute = Reservoir(r).RoutesID % loop on all routes including Rr
        i = Route(iroute).OriginID;
        j = Route(iroute).DestinationID;
        k = Route(iroute).PathID;
        
        Temp_path = Route(iroute).Path;
        if r == i && r == j % internal trip
            Temp_prevresID = i;
            Temp_nextresID = j;
        elseif r == i && r ~= j % if Rr is the route origin
            Temp_prevresID = 0;
            Temp_nextresID = Temp_path(Reservoir(r).RoutesPathIndex(Temp_index)+1); % next reservoir in the path
        elseif r == j && r ~= i % if Rr is the route destination
            Temp_prevresID = Temp_path(Reservoir(r).RoutesPathIndex(Temp_index)-1); % previous reservoir in the path
            Temp_nextresID = 0;
        else
            Temp_prevresID = Temp_path(Reservoir(r).RoutesPathIndex(Temp_index)-1); % previous reservoir in the path
            Temp_nextresID = Temp_path(Reservoir(r).RoutesPathIndex(Temp_index)+1); % next reservoir in the path
        end
        Temp_index = Temp_index + 1;
        
        % Flows for adjacent reservoirs
        for iadj = Reservoir(r).AdjacentRes % loop on all adjacent reservoirs
            if iadj == Temp_prevresID && iadj == Temp_nextresID
                Temp_qin = Reservoir(r).InflowPerODPerPath(i,j,k,:);
                Reservoir(r).InflowPerResPerDest(iadj,j,:) = Reservoir(r).InflowPerResPerDest(iadj,j,:) + Temp_qin(1,1,:);
                Temp_qout = Reservoir(r).OutflowPerODPerPath(i,j,k,:);
                Reservoir(r).OutflowPerResPerDest(iadj,j,:) = Reservoir(r).OutflowPerResPerDest(iadj,j,:) + Temp_qout(1,1,:);
            elseif iadj == Temp_prevresID % if the adjacent reservoir is included in the route
                Temp_qin = Reservoir(r).InflowPerODPerPath(i,j,k,:);
                Reservoir(r).InflowPerResPerDest(iadj,j,:) = Reservoir(r).InflowPerResPerDest(iadj,j,:) + Temp_qin(1,1,:);
            elseif iadj == Temp_nextresID
                Temp_qout = Reservoir(r).OutflowPerODPerPath(i,j,k,:);
                Reservoir(r).OutflowPerResPerDest(iadj,j,:) = Reservoir(r).OutflowPerResPerDest(iadj,j,:) + Temp_qout(1,1,:);
            end
        end
    end
end


for i = 1:(NumTimes-1)
    for r = 1:NumRes
        Reservoir(r).NinPerResPerDest(:,:,i+1) = sum(Reservoir(r).InflowPerResPerDest(:,:,1:i),3)*TimeStep;
        Reservoir(r).NinInternalPerDest(:,i+1) = sum(Reservoir(r).InflowInternalPerDest(:,1:i),2)*TimeStep;
        Reservoir(r).NinTransferPerDest(:,i+1) = sum(Reservoir(r).InflowTransferPerDest(:,1:i),2)*TimeStep;
        Reservoir(r).Nin(i+1) = sum(Reservoir(r).Inflow(1:i))*TimeStep;
        
        Reservoir(r).NoutPerResPerDest(:,:,i+1) = sum(Reservoir(r).OutflowPerResPerDest(:,:,1:i),3)*TimeStep;
        Reservoir(r).NoutPerDest(:,i+1) = sum(Reservoir(r).OutflowPerDest(:,1:i),2)*TimeStep;
        Reservoir(r).Nout(i+1) = sum(Reservoir(r).Outflow(1:i))*TimeStep;
    end
end

% Average trip lengths
for r = 1:NumRes
    Reservoir(r).AvgTripLength = Temp_sumavgtriplength(r,:)./Temp_countavgtriplength(r,:);
end



clear Temp_*

