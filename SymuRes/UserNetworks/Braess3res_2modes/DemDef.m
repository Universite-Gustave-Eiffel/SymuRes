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


%% Scenario 1: demand peak with 3 car and 1 PT routes
%--------------------------------------------------------------------------
% Test scenario from Mahendra Paipuri

if strcmp(Simulation.Name(1:4),'SC11')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.7; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC1
end

if strcmp(Simulation.Name(1:4),'SC12')
    Simulation.MergeModel = 'endogenous'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.7; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC1
end

if strcmp(Simulation.Name(1:4),'SC13')
    Simulation.MergeModel = 'demfifo'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 1.0; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC1
end

if strcmp(Simulation.Name(1:4),'SC14')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'queuedyn'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.7; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC1
end


%% Scenario 2: congestion with 2 car and 1 PT routes, bottleneck at R3 exit
%--------------------------------------------------------------------------
% Test inspired by the scenario SC1 from the Braess network with 1 mode

if strcmp(Simulation.Name(1:4),'SC21')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC2
end


%% Scenario 3: congestion with 3 car routes, bottleneck at R3 exit
%--------------------------------------------------------------------------
% Test inspired by the scenario SC1 from the Braess network with 1 mode

if strcmp(Simulation.Name(1:4),'SC31')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC3
end


rmpath(['UserNetworks/' Simulation.Network '/scenarios/'])

