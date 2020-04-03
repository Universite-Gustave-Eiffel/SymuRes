%% DEMAND DEFINITION
%--------------------------------------------------------------------------

% All macro OD pairs
%--------------------------------------------------------------------------
for od = 1:NumODmacro
    i = 1;
    ODmacro(od).Demand(i).Purpose = 'cartrip';
    ODmacro(od).Demand(i).Time = 0; % [s]
    ODmacro(od).Demand(i).Data = 0; % [veh/s]
end


%% Scenario 1: predefined routes
%--------------------------------------------------------------------------

% Assignment.Periods = [0 Simulation.Duration];
% Assignment.PredefRoute = 1;
% Assignment.Convergence = 0;
% 
% % [R1 R2 R4]
% %--------------------------------------------------------------------------
% iroute = 1;
% for od = 1:NumODmacro
%     if ODmacro(od).NodeOriginID == 1 && ODmacro(od).NodeDestinationID == 3
%         od0 = od;
%     end
% end
% Route(iroute).ODmacroID = od0;
% Route(iroute).ResPath = ODmacro(od0).PossibleRoute(1).ResPath;
% Route(iroute).NodePath = ODmacro(od0).PossibleRoute(1).NodePath;
% Route(iroute).TripLengths = ODmacro(od0).PossibleRoute(1).TripLengths;
% j = 1;
% Route(iroute).Demand0(j).Purpose = 'cartrip';
% Td1 = 500;
% Td11 = 1000;
% Td12 = 2000;
% Td2 = 2500;
% q0 = 0.35;
% q1 = 1.1;
% q2 = 0.35;
% % qinD = @(t_) q0*bumpfct(t_,Td1,Td2,q1/q0,q2/q0);
% qinD = @(t_) q0*trayfct(t_,Td1,Td11,Td12,Td2,q1/q0,q2/q0);
% Route(iroute).Demand0(j).Time = [0 Td1:60:(Td2+60)]; % [s]
% Route(iroute).Demand0(j).Data = [q0 qinD(Td1:60:Td2) q2]; % [veh/s]
% 
% % [R1 R2 R3 R4]
% %--------------------------------------------------------------------------
% iroute = 2;
% for od = 1:NumODmacro
%     if ODmacro(od).NodeOriginID == 1 && ODmacro(od).NodeDestinationID == 3
%         od0 = od;
%     end
% end
% Route(iroute).ODmacroID = od0;
% Route(iroute).ResPath = ODmacro(od0).PossibleRoute(3).ResPath;
% Route(iroute).NodePath = ODmacro(od0).PossibleRoute(3).NodePath;
% Route(iroute).TripLengths = ODmacro(od0).PossibleRoute(3).TripLengths;
% j = 1;
% Route(iroute).Demand0(j).Purpose = 'cartrip';
% Td1 = 500;
% Td11 = 1000;
% Td12 = 1500;
% Td2 = 2000;
% q0 = 0.2;
% q1 = 0.9; %1.2;
% q2 = 0.2;
% % qinD = @(t_) q0*bumpfct(t_,Td1,Td2,q1/q0,q2/q0);
% qinD = @(t_) q0*trayfct(t_,Td1,Td11,Td12,Td2,q1/q0,q2/q0);
% Route(iroute).Demand0(j).Time = [0 Td1:60:(Td2+60)]; % [s]
% Route(iroute).Demand0(j).Data = [q0 qinD(Td1:60:Td2) q2]; % [veh/s]
% 
% % [R1 R3 R4]
% %--------------------------------------------------------------------------
% iroute = 3;
% for od = 1:NumODmacro
%     if ODmacro(od).NodeOriginID == 1 && ODmacro(od).NodeDestinationID == 3
%         od0 = od;
%     end
% end
% Route(iroute).ODmacroID = od0;
% Route(iroute).ResPath = ODmacro(od0).PossibleRoute(2).ResPath;
% Route(iroute).NodePath = ODmacro(od0).PossibleRoute(2).NodePath;
% Route(iroute).TripLengths = ODmacro(od0).PossibleRoute(2).TripLengths;
% j = 1;
% Route(iroute).Demand0(j).Purpose = 'cartrip';
% Td1 = 500;
% Td11 = 1000;
% Td12 = 2000;
% Td2 = 2500;
% q0 = 0.1;
% q1 = 0.2;
% q2 = 0.1;
% % qinD = @(t_) q0*bumpfct(t_,Td1,Td2,q1/q0,q2/q0);
% qinD = @(t_) q0*trayfct(t_,Td1,Td11,Td12,Td2,q1/q0,q2/q0);
% Route(iroute).Demand0(j).Time = [0 Td1:60:(Td2+60)]; % [s]
% Route(iroute).Demand0(j).Data = [q0 qinD(Td1:60:Td2) q2]; % [veh/s]
% 
% MacroNode(9).Capacity.Time = [0 1000];
% MacroNode(9).Capacity.Data = [10 0.4];


%% Scenario 2: one assignment period
%--------------------------------------------------------------------------

% Assignment.Periods = [0 Simulation.Duration];
% Assignment.Convergence = 1;
% Assignment.PredefRoute = 0;
% Assignment.model = 1; % DUE
% Assignment.Behavior = 1; % rational
% 
% % R1 > R4
% %--------------------------------------------------------------------------
% for od = 1:NumODmacro
%     if ODmacro(od).NodeOriginID == 1 && ODmacro(od).NodeDestinationID == 3
%         od0 = od;
%     end
% end
% i = 1;
% ODmacro(od0).Demand(i).Purpose = 'cartrip';
% Td1 = 1000;
% Td11 = 1500;
% Td12 = 2500;
% Td2 = 3000;
% q0 = 0.5;
% q1 = 1.9; % 1.3;
% q2 = 0.5;
% % qinD = @(t_) q0*bumpfct(t_,Td1,Td2,q1/q0,q2/q0);
% qinD = @(t_) q0*trayfct(t_,Td1,Td11,Td12,Td2,q1/q0,q2/q0);
% ODmacro(od0).Demand(i).Time = [0 Td1:60:(Td2+60)]; % [s]
% ODmacro(od0).Demand(i).Data = [q0 qinD(Td1:60:Td2) q2]; % [veh/s]


%% Scenario 3: several assignment periods
%--------------------------------------------------------------------------

Assignment.Periods = [0 1000 2000 3000 4000 Simulation.Duration];
Assignment.Convergence = 1;
Assignment.PredefRoute = 0;
Assignment.model = 1; % DUE
Assignment.Behavior = 1; % rational

% R1 > R4
%--------------------------------------------------------------------------
for od = 1:NumODmacro
    if ODmacro(od).NodeOriginID == 1 && ODmacro(od).NodeDestinationID == 3
        od0 = od;
    end
end
i = 1;
ODmacro(od0).Demand(i).Purpose = 'cartrip';
Td1 = 1000;
Td11 = 1500;
Td12 = 2500;
Td2 = 3000;
q0 = 0.5;
q1 = 1.9; % 1.3;
q2 = 0.5;
% qinD = @(t_) q0*bumpfct(t_,Td1,Td2,q1/q0,q2/q0);
qinD = @(t_) q0*trayfct(t_,Td1,Td11,Td12,Td2,q1/q0,q2/q0);
ODmacro(od0).Demand(i).Time = [0 Td1:60:(Td2+60)]; % [s]
ODmacro(od0).Demand(i).Data = [q0 qinD(Td1:60:Td2) q2]; % [veh/s]


