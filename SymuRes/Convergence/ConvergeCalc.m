%% CONVERGECALC: Convergence calculation
%--------------------------------------------------------------------------
% This routine calculates the convergence criteria. This includes the Gap
% function and the number of violations.
%
% Feb 2019 - Guilhem Mariotte & Sergio F. A. Batista

SimulationDuration = Simulation.Duration;
TimeStep = Simulation.TimeStep;
CurrentTime = Simulation.Time;
Temp_iter = Assignment.CurIteration;

% Assignment period time window
Temp_StartTimeID = floor(Assignment.CurrentTime/TimeStep) + 1;
Temp_EndTimeID = min([floor(Assignment.Periods(Assignment.CurrentPeriodID+1)/TimeStep) NumTimes-1]);

% Initialize the variables:
Temp_Gap = 0;
Temp_NumViolations = 0;

% Loop over all routes to calculate the route utilities:
Temp_Utility = zeros(1,NumRoutes);
for iroute = 1:NumRoutes
    i_p = 1;
    for r = Route(iroute).ResPath
        Temp_meanV = mean(Reservoir(r).MeanSpeed(Temp_StartTimeID:Temp_EndTimeID));
        if Temp_meanV > 0
            Temp_Utility(iroute) = Temp_Utility(iroute) + Route(iroute).TripLengths(i_p)/Temp_meanV;
        else
            Temp_Utility(iroute) = Temp_Utility(iroute) + Route(iroute).TripLengths(i_p)/Reservoir(r).FreeflowSpeed;
        end
        i_p = i_p + 1;
    end
end

% Calculate the Gap:
% Gap for DUE and SUE:
if Assignment.Behavior == 1
    for od = 1:length(ODmacro) % loop on all OD
        Temp_demdata = 0;
        for k = 1:length(ODmacro(od).Demand)
            Temp_demdata = Temp_demdata + sum(ODmacro(od).Demand(k).Data);
        end
        if Temp_demdata > 0 % if the OD (i,j) is used at one time at least
            Temp_RouteIDs = ODmacro(od).RoutesID;
            Temp_coeffs = Assignment.NewAssignCoeff(Temp_RouteIDs);
            Temp_Utils = Temp_Utility(Temp_RouteIDs);
            Temp_minUtil = min(Temp_Utils);
            %Temp_Gap = Temp_Gap + sum(Temp_coeffs.*ODmacro(od).NVehicules.*(Temp_Utils - Temp_minUtil))./(ODmacro(od).NVehicules.*Temp_minUtil);
            Temp_Gap = Temp_Gap + sum(Temp_coeffs.*(Temp_Utils - Temp_minUtil))./(Temp_minUtil);
        end
    end
    
    % Gap for the BR-DUE and BR-SUE:
elseif Assignment.Behavior == 2
    
    for od = 1:length(ODmacro) % loop on all OD
        Temp_demdata = 0;
        for k = 1:length(ODmacro(od).Demand)
            Temp_demdata = Temp_demdata + sum(ODmacro(od).Demand(k).Data);
        end
        if Temp_demdata > 0 % if the OD (i,j) is used at one time at least
            Temp_RouteIDs = ODmacro(od).RoutesID;
            Temp_coeffs = Assignment.NewAssignCoeff(Temp_RouteIDs);
            Temp_Gap = Temp_Gap + sum(Temp_coeffs.*ODmacro(od).NVehicules.*max([Temp_Utils-Assignment.AL(od) 0]))./(ODmacro(od).NVehicules.*Assignment.AL(od));
        end
    end
    
end

% Calculate the number of violations:
for k = 1:NumRoutes
    if abs(Assignment.NewAssignCoeff(k) - Assignment.OldAssignCoeff(k)) > Assignment.NumViolationsThreshold
        Temp_NumViolations = Temp_NumViolations + 1;
    end
end

% Check the convergence tests:
if Assignment.Convergence == 1 && Assignment.PredefRoute == 0
    if Assignment.ConvergenceIndicatorID == 1
        % The user decides to use the Gap criterion only.
        if Temp_Gap <= Assignment.MinimumGap || Temp_iter >= Assignment.MaxIteration
            Assignment.HasConverged = 1;
        end
    elseif Assignment.ConvergenceIndicatorID == 2
        % The user decides to use the Number of violations criterion only.
        if Temp_NumViolations/NumRoutes <= Assignment.NumViolationsTolerance || Temp_iter >= Assignment.MaxIteration
            Assignment.HasConverged = 1;
        end
    elseif Assignment.ConvergenceIndicatorID == 3
        % The user decides to use the Gap and the Number of violations criteria together.
        if (Temp_Gap <= Assignment.MinimumGap && Temp_NumViolations/NumRoutes <= Assignment.NumViolationsTolerance) || Temp_iter >= Assignment.MaxIteration
            Assignment.HasConverged = 1;
        end
    end
else
    Assignment.HasConverged = 1;
end

disp(['Gap = ' num2str(Temp_Gap,3)])


%% Update Route total travel time
%--------------------------------------------------------------------------

% Calculate the mean travel time of each route
if Simulation.Solver == 1
    
    for iroute = 1:NumRoutes % loop on all routes
        Temp_TT = mean(Route(iroute).TravelTime(Temp_StartTimeID:Temp_EndTimeID));
        if Temp_TT == 0
            Route(iroute).TotalTime = Route(iroute).FreeFlowTravelTime;
        elseif Temp_TT > 0
            Route(iroute).TotalTime = Temp_TT;
        end
    end
    for r = 1:NumRes
        Temp_nin = sum(Reservoir(r).NinPerRoute,1);
        Temp_nout = sum(Reservoir(r).NoutPerRoute,1);
        Reservoir(r).TotalTravelTime = ExperiencedTravelTime(Simulation.Time,Temp_nin,Temp_nout);
        Reservoir(r).TotalTravelTime(1) = 0; % TT not defined at t = 0 if acc ini = 0
    end
    
elseif Simulation.Solver == 2
    
    % To recover the distribution of travel times per reservoir. Then, loop
    % over all the vehicles, since on the trip-based model, we know the
    % entry and exit times of each vehicle on each reservoir.
    for iveh = 1:NumVeh
        % Recovering the travel time of the i-th vehicle on each reservoir.
        Temp_vehTT = Vehicle(iveh).ExitTimes - Vehicle(iveh).EntryTimes;
        % For the path of the i-th vehicle, save the reservoir travel time
        % on the corresponding Reservoir.TotaltravelTime structure.
        i_p = 1;
        for r = Route(Vehicle(iveh).RouteID).ResPath
            Reservoir(r).TotalTravelTime(iveh) = Temp_vehTT(i_p);
            i_p = i_p + 1;
        end
    end
    % Loop over all the reservoir to remove the zeros from the Reservoir.TotaltravelTime structure.
    for r = 1:NumRes
        Reservoir(r).TotalTravelTime = Reservoir(r).TotalTravelTime(Reservoir(r).TotalTravelTime~=0);
        Reservoir(r).TotalTravelTime = Reservoir(r).TotalTravelTime(Reservoir(r).TotalTravelTime~=Inf);
        Reservoir(r).TotalTravelTime(find(isnan(Reservoir(r).TotalTravelTime))) = [];
    end
    for r = 1:NumRes
        if isempty(Reservoir(r).TotalTravelTime)==1
            Reservoir(r).TotalTravelTime = 0;
        end
    end
    % Calculating the mean travel times of each macro-route.
    % Here, we assume that the utilities are additive.
    for iroute = 1:NumRoutes % loop on all routes
        Temp_routeTT = 0;
        for r = Route(iroute).ResPath
            Temp_routeTT = Temp_routeTT + mean(Reservoir(r).TotalTravelTime);
        end
        % Travel time is converted from seconds to hours because of the
        % exponentials in the Logit and Mixed Logit.
        Route(iroute).TotalTime = Temp_routeTT;
    end
    
end


%% Save the simulation results
%--------------------------------------------------------------------------

% Snapshot of the simulation results for the next assignment period
if Assignment.HasConverged == 1
    
    for iroute = 1:NumRoutes
        Route(iroute).ListAssignCoeff = [Route(iroute).ListAssignCoeff Route(iroute).AssignCoeff];
    end
    
    Snapshot.Reservoir = Reservoir;
    Snapshot.Route = Route;
    
    if Simulation.Solver == 2
        Global.CurrentTimeID = Global.NumEvents + 1;
        Snapshot.Global = Global;
        Snapshot.Vehicle = Vehicle;
        Snapshot.NextEvent = NextEvent;
    end
    
    if Assignment.Periods(Assignment.CurrentPeriodID+1) >= Simulation.Duration
        % End of the simulation
        
        % Track eventual gridlocks
        Temp_listres = [];
        for r = 1:NumRes
            if Simulation.Solver == 1
                Temp_acctest = sum(Reservoir(r).Acc >= Reservoir(r).MaxAcc);
            else
                Temp_acctest = 0;
            end
            Temp_speedtest = sum(Reservoir(r).MeanSpeed < 0);
            if Temp_acctest > 0 || Temp_speedtest > 0
                Temp_listres = [Temp_listres r];
            end
        end
        if ~isempty(Temp_listres)
            warning(['The following reservoir(s) encountered gridlock: ' int2str(Temp_listres)])
        end
        
        if Simulation.OutfileDebug == 0
            % Only keep the essential information in the Reservoir and Route structures
            if Simulation.Solver == 1 % accumulation-based solver
                for r = 1:NumRes
                    Out.Reservoir(r).MeanSpeed = Reservoir(r).MeanSpeed;
                    Out.Reservoir(r).AccPerRoute = Reservoir(r).AccPerRoute;
                    Out.Reservoir(r).InflowPerRoute = Reservoir(r).InflowPerRoute;
                    Out.Reservoir(r).OutflowPerRoute = Reservoir(r).OutflowPerRoute;
                    if isfield(Reservoir,'AccCircuPerRoute')
                        Out.Reservoir(r).AccCircuPerRoute = Reservoir(r).AccCircuPerRoute;
                        Out.Reservoir(r).AccQueuePerRoute = Reservoir(r).AccQueuePerRoute;
                        Out.Reservoir(r).OutflowCircuPerRoute = Reservoir(r).OutflowCircuPerRoute;
                    end
                end
                for iroute = 1:NumRoutes
                    Out.Route(iroute).AssignCoeff = Route(iroute).AssignCoeff;
                    Out.Route(iroute).ListAssignCoeff = Route(iroute).ListAssignCoeff;
                    Out.Route(iroute).NVehicles = Route(iroute).NVehicles;
                    Out.Route(iroute).TravelTime = Route(iroute).TravelTime;
                end
            elseif Simulation.Solver == 2 % trip-based solver
                for iroute = 1:NumRoutes
                    Out.Route(iroute).AssignCoeff = Route(iroute).AssignCoeff;
                    Out.Route(iroute).ListAssignCoeff = Route(iroute).ListAssignCoeff;
                    Out.Route(iroute).NVehicles = Route(iroute).NVehicles;
                    Out.Route(iroute).TravelTime = Route(iroute).TravelTime;
                end
            end
            
            Reservoir = Out.Reservoir;
            Route = Out.Route;
        end
    end
end

