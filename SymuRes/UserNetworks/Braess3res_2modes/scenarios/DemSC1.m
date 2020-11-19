%% Scenario 1: demand peak with 3 car and 1 PT routes
%--------------------------------------------------------------------------

Simulation.Duration = 5000; % Simulation duration [s]

Assignment.Periods = [0 Simulation.Duration];
Assignment.Convergence = 0;
Assignment.PredefRoute = 1;

% [R1 R2]
%--------------------------------------------------------------------------
iroute = 1;
for od = 1:NumODmacro
    if ODmacro(od).NodeOriginID == 1 && ODmacro(od).NodeDestinationID == 3
        od0 = od;
    end
end
Route(iroute).ODmacroID = od0;
Route(iroute).ResPath = ODmacro(od0).PossibleRoute(1).ResPath;
Route(iroute).NodePath = ODmacro(od0).PossibleRoute(1).NodePath;
Route(iroute).TripLengths = ODmacro(od0).PossibleRoute(1).TripLengths;
j = 1;
Route(iroute).Demand0(j).Purpose = 'cartrip';
Td1 = 500;
Td2 = 2500;
q0 = 0.05;
q1 = 1.2;
q2 = 0.05;
qinD = @(t_) q0*((t_ <= Td1).*(t_ > 0)) + q1*((Td1 < t_).*(t_ <= Td2)) + q0*(Td2 < t_);
Route(iroute).Demand0(j).Time = [0 Td1:60:(Td2+60)]; % [s]
Route(iroute).Demand0(j).Data = [q0 qinD(Td1:60:Td2) q2]; % [veh/s]

% [R1 R2]
%--------------------------------------------------------------------------
iroute = 2;
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
Td2 = 2500;
q0 = 0.03;
q1 = 0.7; %1.2;
q2 = 0.03;
qinD = @(t_) q0*((t_ <= Td1).*(t_ > 0)) + q1*((Td1 < t_).*(t_ <= Td2)) + q0*(Td2 < t_);
Route(iroute).Demand0(j).Time = [0 Td1:60:(Td2+60)]; % [s]
Route(iroute).Demand0(j).Data = [q0 qinD(Td1:60:Td2) q2]; % [veh/s]

% [R1 R2]
%--------------------------------------------------------------------------
iroute = 3;
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
Td2 = 2500;
q0 = 0.02;
q1 = 0.7;
q2 = 0.02;
qinD = @(t_) q0*((t_ <= Td1).*(t_ > 0)) + q1*((Td1 < t_).*(t_ <= Td2)) + q0*(Td2 < t_);
Route(iroute).Demand0(j).Time = [0 Td1:60:(Td2+60)]; % [s]
Route(iroute).Demand0(j).Data = [q0 qinD(Td1:60:Td2) q2]; % [veh/s]

% [R1 R2]
%--------------------------------------------------------------------------
iroute = 4;
for od = 1:NumODmacro
    if ODmacro(od).NodeOriginID == 1 && ODmacro(od).NodeDestinationID == 3
        od0 = od;
    end
end
Route(iroute).ODmacroID = od0;
Route(iroute).ResPath = ODmacro(od0).PossibleRoute(4).ResPath;
Route(iroute).NodePath = ODmacro(od0).PossibleRoute(4).NodePath;
Route(iroute).TripLengths = ODmacro(od0).PossibleRoute(4).TripLengths;
j = 1;
Route(iroute).Demand0(j).Purpose = 'publictransport';
Td1 = 500;
Td2 = 2500;
q0 = 0.01;
q1 = 0.22; %1.2;
q2 = 0.01;
qinD = @(t_) q0*((t_ <= Td1).*(t_ > 0)) + q1*((Td1 < t_).*(t_ <= Td2)) + q0*(Td2 < t_);
Route(iroute).Demand0(j).Time = [0 Td1:60:(Td2+60)]; % [s]
Route(iroute).Demand0(j).Data = [q0 qinD(Td1:60:Td2) q2]; % [veh/s]

