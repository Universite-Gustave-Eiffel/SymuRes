%% Scenario 3: two trips, demand surge
%--------------------------------------------------------------------------
% Fig. 8 (section 4.1) of Mariotte & Leclercq (Part B 2019)

Assignment.Periods = [0 Simulation.Duration];
Assignment.PredefRoute = 1;
Assignment.Convergence = 0;

% Entry supply function
i = 1;
Reservoir(i).EntryfctParam = [Reservoir(i).MaxAcc Reservoir(i).CritAcc Reservoir(i).MaxProd ...
    0.5*Reservoir(i).CritAcc 1.5*Reservoir(i).CritAcc 1.4*Reservoir(i).MaxProd];

% Entry demand, route 1
%--------------------------------------------------------------------------
iroute = 1;
for od = 1:NumODmacro
    if ODmacro(od).NodeOriginID == 1 && ODmacro(od).NodeDestinationID == 3
        od0 = od;
    end
end
Route(iroute).ODmacroID = od0;
Route(iroute).ResPath = ODmacro(od0).PossibleRoute(2).ResPath;
Route(iroute).NodePath = ODmacro(od0).PossibleRoute(2).NodePath;
Route(iroute).TripLengths = ODmacro(od0).PossibleRoute(2).TripLengths;
j = 1;
Route(iroute).Demand0(j).Purpose = 'cartrip';
Td1 = 1500;
q0 = 0.3;
q1 = 0.8;
Route(iroute).Demand0(j).Time = [0 Td1]; % [s]
Route(iroute).Demand0(j).Data = [q0 q1]; % [veh/s]

% Entry demand, route 2
%--------------------------------------------------------------------------
iroute = 2;
for od = 1:NumODmacro
    if ODmacro(od).NodeOriginID == 2 && ODmacro(od).NodeDestinationID == 4
        od0 = od;
    end
end
Route(iroute).ODmacroID = od0;
Route(iroute).ResPath = ODmacro(od0).PossibleRoute(1).ResPath;
Route(iroute).NodePath = ODmacro(od0).PossibleRoute(1).NodePath;
Route(iroute).TripLengths = ODmacro(od0).PossibleRoute(1).TripLengths;
j = 1;
Route(iroute).Demand0(j).Purpose = 'cartrip';
Td1 = 1500;
q0 = 0.3;
q1 = 1.1;
Route(iroute).Demand0(j).Time = [0 Td1]; % [s]
Route(iroute).Demand0(j).Data = [q0 q1]; % [veh/s]

% Exit supply, route 1
%--------------------------------------------------------------------------
i = 3;
MacroNode(i).Capacity.Time = 0; % [s]
MacroNode(i).Capacity.Data = 0.6; % [veh/s]

% Exit supply, route 2
%--------------------------------------------------------------------------
i = 4;
MacroNode(i).Capacity.Time = 0; % [s]
MacroNode(i).Capacity.Data = 0.7; % [veh/s]


