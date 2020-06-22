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

addpath(['UserNetworks/' Simulation.Network '/scenarios/'])


%% Scenario 1: uniform loading
%--------------------------------------------------------------------------

if strcmp(Simulation.Name(1:4),'SC11')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC1
end


%% Scenario 2: loading with a gravity model
%--------------------------------------------------------------------------

if strcmp(Simulation.Name(1:4),'SC21')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC2
end



rmpath(['UserNetworks/' Simulation.Network '/scenarios/'])

