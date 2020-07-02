%% ASSIGNCALC: Assignment calculation
%--------------------------------------------------------------------------
% Assign OD demand flows to the routes (path flow distribution):
% Calculate the new possible assignment or path flow coefficients for each route
% with a choice model, and update the effective coefficient with the MSA formula.
%
% June 2020 - Guilhem Mariotte

SimulationDuration = Simulation.Duration;
TimeStep = Simulation.TimeStep;

% Print route info
fprintf('%s \n','Travel time per route:')
for iroute = 1:NumRoutes
    fprintf('%3.1f \t',Route(iroute).TotalTime)
end

% Load from previous assignment period
if Assignment.CurrentPeriodID > 1
    Route = Snapshot.Route;
end

Assignment.OldAssignCoeff = Assignment.NewAssignCoeff;

% Calculation of the MSWA parameter at the current iteration
iter = Assignment.CurIteration;
Temp_betaMSWA = Assignment.CurIteration^Assignment.WeightMSWA;
Assignment.GammaMSWA = Assignment.GammaMSWA + Temp_betaMSWA;
Temp_alphaMSWA = Temp_betaMSWA/Assignment.GammaMSWA;

% Initializing the AL in the structure assignment for BR-DUE or BR-SUE.
if Assignment.model == 3 || Assignment.model == 4
    Assignment.AL = zeros(1,NumODmacro);
end

if Assignment.PredefRoute == 0
    % If the routes are calculated by the assignment procedure
    %---------------------------------------------------------
    
    for od = 1:NumODmacro % loop on all OD
        % List of routes for the OD (o,d)
        Temp_RouteIDs = ODmacro(od).RoutesID;
        if ~isempty(Temp_RouteIDs)
            % Apply a choice model to calculate the path flow distribution (assignment coefficients)
            if Assignment.model < 100
                % Choice model
                ODflowdistrib = ChoiceModel(od,Temp_RouteIDs,Reservoir,Route,Assignment,Simulation);
            else
                % Empirical model
                ODflowdistrib = EmpiricalModel(od,Temp_RouteIDs,Reservoir,Route,Assignment,Simulation);
            end
            Temp_coeffs = ODflowdistrib(Temp_RouteIDs);
            Temp_oldcoeffs = Assignment.OldAssignCoeff(Temp_RouteIDs);
            
            % MSA step to update the assignment coefficients
            % Classical MSA
            Assignment.NewAssignCoeff(Temp_RouteIDs) = (1 - 1/Assignment.CurIteration).*Temp_oldcoeffs + 1/Assignment.CurIteration.*Temp_coeffs;
            % Weighted MSA (MSWA)
            Assignment.NewAssignCoeff(Temp_RouteIDs) = Temp_oldcoeffs + Temp_alphaMSWA.*(Temp_coeffs - Temp_oldcoeffs);
            
            for iroute = Temp_RouteIDs % loop on all routes for the OD (i,j)
                % Assignment coefficient
                Route(iroute).AssignCoeff = Assignment.NewAssignCoeff(iroute);
                % Route demand
                for j = 1:length(ODmacro(od).Demand)
                    if strcmp(ODmacro(od).Demand(j).Purpose,'cartrip')
                        Route(iroute).Demand = Route(iroute).AssignCoeff.*ODmacro(od).Demand(j).Data2; % row vector
                    end
                end
                % Append the number of vehicles assigned to each route
                Route(iroute).NVehicles = ODmacro(od).NVehicules*Route(iroute).AssignCoeff;
            end
        end
    end
    
else
    % If the routes with their demand are given prior to the simulation
    %------------------------------------------------------------------
    
    for od = 1:NumODmacro
        Temp_demtot = 0;
        for iroute = ODmacro(od).RoutesID
            Temp_demtot = Temp_demtot + sum(Route(iroute).Demand);
        end
        for iroute = ODmacro(od).RoutesID
            Route(iroute).AssignCoeff = sum(Route(iroute).Demand)./Temp_demtot;
        end
    end
    for iroute = 1:NumRoutes
        od = Route(iroute).ODmacroID;
        Route(iroute).NVehicles = ODmacro(od).NVehicules;
        Assignment.NewAssignCoeff(iroute) = Route(iroute).AssignCoeff;
    end
end

% Set small assignment coefficients to zero to avoid numerical issues
Temp_coeff_threshold = 0.005;
for od = 1:NumODmacro
    Temp_coefftot = 0;
    for iroute = ODmacro(od).RoutesID
        if Route(iroute).AssignCoeff < Temp_coeff_threshold
            Route(iroute).AssignCoeff = 0;
        end
        Temp_coefftot = Temp_coefftot + Route(iroute).AssignCoeff;
    end
    for iroute = ODmacro(od).RoutesID
        Route(iroute).AssignCoeff = Route(iroute).AssignCoeff./Temp_coefftot;
        Assignment.NewAssignCoeff(iroute) = Route(iroute).AssignCoeff;
    end
end


% Print route info
fprintf('\n%s \n','Assignment coefficient per route:')
for iroute = 1:NumRoutes
    fprintf('%1.3f \t',Route(iroute).AssignCoeff)
end
fprintf('%s\n',' ')


%% Initial average trip length estimation
%--------------------------------------------------------------------------

if Assignment.CurrentPeriodID == 1
    for r = 1:NumRes
        Reservoir(r).AvgTripLength = zeros(NumModes,NumTimes);
        Temp_demtot = zeros(1,NumModes);
        Temp_demLtrip = zeros(1,NumModes);
        i_r = 1;
        for iroute = Reservoir(r).RoutesID
            od = Route(iroute).ODmacroID;
            i_m = Route(iroute).ModeID;
            Temp_dem = Route(iroute).Demand(1);
            Temp_Ltrip = Reservoir(r).TripLengthPerRoute(i_r);
            
            Temp_demtot(i_m) = Temp_demtot(i_m) + Temp_dem;
            Temp_demLtrip(i_m) = Temp_demLtrip(i_m) + Temp_Ltrip*Temp_dem;
            i_r = i_r + 1;
        end
        for i_m = 1:NumModes
            if Temp_demtot(i_m) > 0
                Reservoir(r).AvgTripLength(i_m,1) = Temp_demLtrip(i_m)/Temp_demtot(i_m);
            else
                Reservoir(r).AvgTripLength(i_m,1) = mean(Reservoir(r).TripLengthPerRoute);
            end
        end
    end
end

