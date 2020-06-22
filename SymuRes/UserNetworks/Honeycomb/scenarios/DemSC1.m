%% Scenario 1: uniform loading
%--------------------------------------------------------------------------

Simulation.Duration = 5000; % Simulation duration [s]

Assignment.Periods = [0 Simulation.Duration];
Assignment.PredefRoute = 0;
Assignment.Convergence = 0;

% Peak profile
Td1 = 1000;
Td11 = 1500;
Td12 = 2500;
Td2 = 3000;
q0 = 0.01;
q1 = 0.04;
q2 = 0.01;
qinD = @(t_) q0*trayfct(t_,Td1,Td11,Td12,Td2,q1/q0,q2/q0);

% Demand for all macro OD
%--------------------------------------------------------------------------
for od = 1:NumODmacro
    i = 1;
    ODmacro(od).Demand(i).Purpose = 'cartrip';
    ODmacro(od).Demand(i).Time = [0 Td1:60:(Td2+60)]; % [s]
    ODmacro(od).Demand(i).Data = [q0 qinD(Td1:60:Td2) q2]; % [veh/s]
end

