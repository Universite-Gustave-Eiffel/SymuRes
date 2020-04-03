%% Scenario 2: one trip, supply reduction
%--------------------------------------------------------------------------
% Fig. 4.a2.b2.c2 (section 2.3) of Mariotte & Leclercq (Part B 2019)

Assignment.Periods = [0 Simulation.Duration];
Assignment.PredefRoute = 1;
Assignment.Convergence = 0;

% Entry supply function
i = 1;
Reservoir(i).EntryfctParam = [Reservoir(i).MaxAcc Reservoir(i).CritAcc Reservoir(i).MaxProd ...
    0.8*Reservoir(i).CritAcc 1*Reservoir(i).CritAcc 1*Reservoir(i).MaxProd];

% Entry demand
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
Route(iroute).Demand0(j).Time = 0; % [s]
Route(iroute).Demand0(j).Data = 0.7; % [veh/s]

% Exit supply
%--------------------------------------------------------------------------
i = 3;
Td1 = 1000;
Td2 = 6000;
q0 = 1.2;
q1 = 0.1;
q2 = 1.2;
qoutS = @(t_) q0*bumpfct(t_,Td1,Td2,q1/q0,q2/q0);
MacroNode(i).Capacity.Time = [0 Td1:20:(Td2+20)]; % [s]
MacroNode(i).Capacity.Data = [q0 qoutS(Td1:20:Td2) q2]; % [veh/s]


