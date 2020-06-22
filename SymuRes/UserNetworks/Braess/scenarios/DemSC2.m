%% Scenario 2: congestion on predefined routes, bottleneck at R4 exit
%--------------------------------------------------------------------------

Simulation.Duration = 8000; % Simulation duration [s]

Assignment.Periods = [0 Simulation.Duration];
Assignment.PredefRoute = 1;
Assignment.Convergence = 0;

% Demand route [R1 R3 R4]
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
Td1 = 500;
Td11 = 1000;
Td12 = 2000;
Td2 = 2500;
q0 = 0.1;
q1 = 0.8;
q2 = 0.1;
qinD = @(t_) q0*trayfct(t_,Td1,Td11,Td12,Td2,q1/q0,q2/q0);
Route(iroute).Demand0(j).Time = [0 Td1:60:(Td2+60)]; % [s]
Route(iroute).Demand0(j).Data = [q0 qinD(Td1:60:Td2) q2]; % [veh/s]

% Demand route [R1 R2 R3 R4]
%--------------------------------------------------------------------------
iroute = 2;
for od = 1:NumODmacro
    if ODmacro(od).NodeOriginID == 1 && ODmacro(od).NodeDestinationID == 3
        od0 = od;
    end
end
Route(iroute).ODmacroID = od0;
Route(iroute).ResPath = ODmacro(od0).PossibleRoute(3).ResPath;
Route(iroute).NodePath = ODmacro(od0).PossibleRoute(3).NodePath;
Route(iroute).TripLengths = ODmacro(od0).PossibleRoute(3).TripLengths;
j = 1;
Route(iroute).Demand0(j).Purpose = 'cartrip';
Td1 = 500;
Td11 = 1000;
Td12 = 1500;
Td2 = 2000;
q0 = 0.2;
q1 = 0.9; %1.2;
q2 = 0.2;
qinD = @(t_) q0*trayfct(t_,Td1,Td11,Td12,Td2,q1/q0,q2/q0);
Route(iroute).Demand0(j).Time = [0 Td1:60:(Td2+60)]; % [s]
Route(iroute).Demand0(j).Data = [q0 qinD(Td1:60:Td2) q2]; % [veh/s]


% Supply exit from R4
%--------------------------------------------------------------------------
MacroNode(3).Type = 'externalexit';
MacroNode(3).Coord = Reservoir(4).Centroid + [1 0.2];
MacroNode(3).Capacity.Time = [0 1000 6000];
MacroNode(3).Capacity.Data = [10 0.4 10];


