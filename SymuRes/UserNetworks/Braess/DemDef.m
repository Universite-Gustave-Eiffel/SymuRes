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


%% Scenario 1: congestion on predefined routes, bottleneck at R3 exit
%--------------------------------------------------------------------------

if strcmp(Simulation.Name(1:4),'SC11')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC1
end

if strcmp(Simulation.Name(1:4),'SC12')
    Simulation.MergeModel = 'endogenous'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC1
end

if strcmp(Simulation.Name(1:4),'SC13')
    Simulation.MergeModel = 'demfifo'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC1
end


%% Scenario 2: congestion on predefined routes, bottleneck at R4 exit
%--------------------------------------------------------------------------

if strcmp(Simulation.Name(1:4),'SC21')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.4; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC2
end


%% Scenario 3: DTA with several assignment periods
%--------------------------------------------------------------------------
% ! Low demand, thus does not work with a simufactor < 0.6 in the trip-based

if strcmp(Simulation.Name(1:4),'SC31')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.7; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC3
end


%% Scenario 4: DTA with one assignment period
%--------------------------------------------------------------------------

if strcmp(Simulation.Name(1:4),'SC41')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.7; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC4
end


%% Scenario 5: DTA with heavy congestion
%--------------------------------------------------------------------------

if strcmp(Simulation.Name(1:4),'SC51')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.7; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC5
end


%% Scenario 6: Predefined assignment coeff with heavy congestion
%--------------------------------------------------------------------------

if strcmp(Simulation.Name(1:4),'SC61')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    Simulation.TripbasedSimuFactor = 0.7; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
    DemSC6
end


rmpath(['UserNetworks/' Simulation.Network '/scenarios/'])

