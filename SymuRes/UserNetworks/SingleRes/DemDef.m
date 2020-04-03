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


%% Scenario 1: one trip, demand peak
%--------------------------------------------------------------------------
% Fig. 4.a1.b1.c1 (section 2.3) of Mariotte & Leclercq (Part B 2019)

if strcmp(Simulation.Name,'SC11')
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    DemSC1
end

if strcmp(Simulation.Name,'SC12')
    Simulation.DivergeModel = 'decrdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    DemSC1
end


%% Scenario 2: one trip, supply reduction
%--------------------------------------------------------------------------
% Fig. 4.a2.b2.c2 (section 2.3) of Mariotte & Leclercq (Part B 2019)

if strcmp(Simulation.Name,'SC21')
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    DemSC2
end

if strcmp(Simulation.Name,'SC22')
    Simulation.DivergeModel = 'decrdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    DemSC2
end


%% Scenario 3: two trips, demand surge
%--------------------------------------------------------------------------
% Fig. 8 (section 4.1) of Mariotte & Leclercq (Part B 2019)

if strcmp(Simulation.Name,'SC31')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'decrdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    DemSC3
end

if strcmp(Simulation.Name,'SC32')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    DemSC3
end

if strcmp(Simulation.Name,'SC33')
    Simulation.MergeModel = 'endogenous'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    DemSC3
end


%% Scenario 4: two trips, demand peak
%--------------------------------------------------------------------------
% Not presented in the paper

if strcmp(Simulation.Name,'SC41')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'decrdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    DemSC4
end

if strcmp(Simulation.Name,'SC42')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    DemSC4
end

if strcmp(Simulation.Name,'SC43')
    Simulation.MergeModel = 'endogenous'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    DemSC4
end


%% Scenario 5: two trips, supply reduction
%--------------------------------------------------------------------------
% Fig. 9 (section 4.2) of Mariotte & Leclercq (Part B 2019)

if strcmp(Simulation.Name,'SC51')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'decrdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    DemSC5
end

if strcmp(Simulation.Name,'SC52')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    DemSC5
end

if strcmp(Simulation.Name,'SC53')
    Simulation.MergeModel = 'endogenous'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
    DemSC5
end


rmpath(['UserNetworks/' Simulation.Network '/scenarios/'])

