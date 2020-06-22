%% Scenario 5: DTA with heavy congestion
%--------------------------------------------------------------------------

Simulation.Duration = 7000; % Simulation duration [s]

Assignment.Periods = [0 2000 4000 6000 Simulation.Duration];
% Assignment.Periods = [0 Simulation.Duration];
Assignment.Convergence = 1;
Assignment.PredefRoute = 0;
Assignment.MaxIteration = 15; % Maximum number of iterations for the convergence loop
Assignment.model = 1; % DUE
Assignment.Behavior = 1; % rational

% Demand OD R1 > R4
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
qinD = @(t_) q0*trayfct(t_,Td1,Td11,Td12,Td2,q1/q0,q2/q0);
ODmacro(od0).Demand(i).Time = [0 Td1:60:(Td2+60)]; % [s]
ODmacro(od0).Demand(i).Data = [q0 qinD(Td1:60:Td2) q2]; % [veh/s]

% Demand R2 (internal trips)
%--------------------------------------------------------------------------
for od = 1:NumODmacro
    if ODmacro(od).NodeOriginID == 10 && ODmacro(od).NodeDestinationID == 11
        od0 = od;
    end
end
i = 1;
ODmacro(od0).Demand(i).Purpose = 'cartrip';
ODmacro(od0).Demand(i).Time = 0; % [s]
ODmacro(od0).Demand(i).Data = 0.8; % [veh/s]

% Supply transfer from R2 to R4
%--------------------------------------------------------------------------
MacroNode(8).Capacity.Time = 0;
MacroNode(8).Capacity.Data = 1;

% Supply transfer from R3 to R4
%--------------------------------------------------------------------------
MacroNode(9).Capacity.Time = 0;
MacroNode(9).Capacity.Data = 0.5;

