%% RESERVOIR DEFINITION
%--------------------------------------------------------------------------

%%% RESERVOIR
% Structure which defines all information related to the reservoirs

% Reservoir.Centroid ; virtual positioning [x y] for plotting purpose
% Reservoir.AdjacentRes ; vector, list of adjacent reservoir IDs
% Reservoir.NetLength ; total network length in the reservoir [m]
% Reservoir.FreeflowSpeed ; free-flow speed in the reservoir [m/s]
% Reservoir.MaxProd ; reservoir maximum (critical) production [veh.m/s]
% Reservoir.MaxAcc ; reservoir maximum (jam) accumulation [veh]
% Reservoir.CritAcc ; reservoir critical accumulation [veh]

% Number of reservoirs
NumRes = 3;

Reservoir = struct('Centroid',cell(1,NumRes));

% MFD function
MFDfct = @parabo3dFD;

% Entry supply function
Entryfct = @parabo3dEntryFD;

% Exit demand function
Exitfct = @parabo3dExitFD;

% Reservoir definition

% R1
%--------------------------------------------------------------------------
i = 1;
Reservoir(i).Centroid = [0 1]; % virtual positioning [x y] for plotting purpose
Reservoir(i).BorderPoints = [-1 1 1 -1 -1; 0 0 2 2 0];
Reservoir(i).AdjacentRes = [2 3]; % adjacent reservoirs
Reservoir(i).FreeflowSpeed = [15 15]; % [m/s]
Reservoir(i).MaxProd = 5000; % 3750; % [veh.m/s]
Reservoir(i).CritAcc(1) = 2*Reservoir(i).MaxProd/Reservoir(i).FreeflowSpeed(1); % for a parabolic MFD
Reservoir(i).MaxAcc = 2*Reservoir(i).CritAcc(1); % for a parabolic MFD
% 3DMFD parameters
uc       = Reservoir(i).FreeflowSpeed(1); % Free flow speed of cars
ub       = Reservoir(i).FreeflowSpeed(2); % Free flow speed of buses
bcc      = -Reservoir(i).MaxProd/Reservoir(i).CritAcc^2; % Marginal Effect of cars on car speed
bbc      = -0.3;     % Marginal Effect of buses on car speed
bcb      = 0.2*bcc;  % Marginal Effect of cars on bus speed
bbb      = 0.2*bbc;  % Marginal Effect of buses on bus speed
Reservoir(i).MFDfctParam = [uc ub bcc bbc bcb bbb 1 1];
Reservoir(i).CritAcc(2) = -uc/(bbc + bcb);
Reservoir(i).EntryfctParam = Reservoir(i).MFDfctParam;
Reservoir(i).ExitfctParam = Reservoir(i).MFDfctParam;


% R2
%--------------------------------------------------------------------------
i = 2;
Reservoir(i).Centroid = [2 2]; % virtual positioning [x y]
Reservoir(i).BorderPoints = [1 3 3 1 1; 1 1 3 3 1];
Reservoir(i).AdjacentRes = [1 3 4]; % adjacent reservoirs
Reservoir(i).FreeflowSpeed = [15 15]; % [m/s]
Reservoir(i).MaxProd = 4500; % 2250; % [veh.m/s]
Reservoir(i).CritAcc(1) = 2*Reservoir(i).MaxProd/Reservoir(i).FreeflowSpeed(1); % for a parabolic MFD
Reservoir(i).MaxAcc = 2*Reservoir(i).CritAcc(1); % for a parabolic MFD
% 3DMFD parameters
uc       = Reservoir(i).FreeflowSpeed(1); % Free flow speed of cars
ub       = Reservoir(i).FreeflowSpeed(2); % Free flow speed of buses
bcc      = -Reservoir(i).MaxProd/Reservoir(i).CritAcc^2; % Marginal Effect of cars on car speed
bbc      = -0.3;     % Marginal Effect of buses on car speed
bcb      = 0.2*bcc;  % Marginal Effect of cars on bus speed
bbb      = 0.2*bbc;  % Marginal Effect of buses on bus speed
Reservoir(i).MFDfctParam = [uc ub bcc bbc bcb bbb 1 1];
Reservoir(i).CritAcc(2) = -uc/(bbc + bcb);
Reservoir(i).EntryfctParam = Reservoir(i).MFDfctParam;
Reservoir(i).ExitfctParam = Reservoir(i).MFDfctParam;

% R3
%--------------------------------------------------------------------------
i = 3;
Reservoir(i).Centroid = [4 1]; % virtual positioning [x y]
Reservoir(i).BorderPoints = [3 5 5 3 3; 0 0 2 2 0];
Reservoir(i).AdjacentRes = [2 3]; % adjacent reservoirs
Reservoir(i).FreeflowSpeed = [15 15]; % [m/s]
Reservoir(i).MaxProd = 3000; % 3750; % [veh.m/s]
Reservoir(i).CritAcc(1) = 2*Reservoir(i).MaxProd/Reservoir(i).FreeflowSpeed(1); % for a parabolic MFD
Reservoir(i).MaxAcc = 2*Reservoir(i).CritAcc(1); % for a parabolic MFD
% 3DMFD parameters
uc       = Reservoir(i).FreeflowSpeed(1); % Free flow speed of cars
ub       = Reservoir(i).FreeflowSpeed(2); % Free flow speed of buses
bcc      = -Reservoir(i).MaxProd/Reservoir(i).CritAcc^2; % Marginal Effect of cars on car speed
bbc      = -0.3;     % Marginal Effect of buses on car speed
bcb      = 0.2*bcc;  % Marginal Effect of cars on bus speed
bbb      = 0.2*bbc;  % Marginal Effect of buses on bus speed
Reservoir(i).MFDfctParam = [uc ub bcc bbc bcb bbb 1 1];
Reservoir(i).CritAcc(2) = -uc/(bbc + bcb);
Reservoir(i).EntryfctParam = Reservoir(i).MFDfctParam;
Reservoir(i).ExitfctParam = Reservoir(i).MFDfctParam;

%%% MACRONODE
% Structure that encompasses all the information about the reservoir
% interfaces. Depending on the simulation purpose, any number of nodes
% can be created in each reservoir. They are of 3 types: internal origin,
% internal destination, or interface border with an adjacent reservoir.

% MacroNode.Type ; 'origin', 'destination' or 'border'
% MacroNode.ResID ; ID of the reservoir, or [ID_res1 ID_res2] in case of a
% border type
% MacroNode.Coord ; [x y] location of the MacroNode element for plotting purpose
% MacroNode.Capacity.Time ; row vector of time when the capacity values change
% MacroNode.Capacity.Data ; row vector of capacity values across time

i = 1;
MacroNode(i).Type = 'externalentry';
MacroNode(i).ResID = 1;
MacroNode(i).Coord = Reservoir(1).Centroid + [0 0.2];
MacroNode(i).Capacity.Time = 0; % [s]
MacroNode(i).Capacity.Data = 100; % [veh/s]

i = 2;
MacroNode(i).Type = 'origin';
MacroNode(i).ResID = 1;
MacroNode(i).Coord = Reservoir(1).Centroid + [0 -0.2];
MacroNode(i).Capacity.Time = 0; % [s]
MacroNode(i).Capacity.Data = 100; % [veh/s]

i = 3;
MacroNode(i).Type = 'externalexit';
MacroNode(i).ResID = 3;
MacroNode(i).Coord = Reservoir(3).Centroid + [0 0.2];
MacroNode(i).Capacity.Time = 0; % [s]
MacroNode(i).Capacity.Data = 100; % [veh/s]

i = 4;
MacroNode(i).Type = 'destination';
MacroNode(i).ResID = 3;
MacroNode(i).Coord = Reservoir(3).Centroid + [0 -0.2];
MacroNode(i).Capacity.Time = 0; % [s]
MacroNode(i).Capacity.Data = 100; % [veh/s]

i = 5;
MacroNode(i).Type = 'border';
MacroNode(i).ResID = [1 2];
MacroNode(i).Coord = mean([Reservoir(1).Centroid; Reservoir(2).Centroid]);
MacroNode(i).Capacity.Time = 0; % [s]
MacroNode(i).Capacity.Data = 100; % [veh/s]

i = 6;
MacroNode(i).Type = 'border';
MacroNode(i).ResID = [2 3];
MacroNode(i).Coord = mean([Reservoir(2).Centroid; Reservoir(3).Centroid]);
MacroNode(i).Capacity.Time = 0; % [s]
MacroNode(i).Capacity.Data = 100; % [veh/s]

i = 7;
MacroNode(i).Type = 'border';
MacroNode(i).ResID = [2 3];
MacroNode(i).Coord = mean([Reservoir(2).Centroid - 0.1; Reservoir(3).Centroid - 0.1]);
MacroNode(i).Capacity.Time = 0; % [s]
MacroNode(i).Capacity.Data = 100; % [veh/s]


NumMacroNodes = i;

% Append to Reservoir structure
for r = 1:NumRes
    Reservoir(r).MacroNodesID = [];
    for i = 1:NumMacroNodes
        if ismember(r,MacroNode(i).ResID)
            Reservoir(r).MacroNodesID = [Reservoir(r).MacroNodesID i];
        end
    end
    Reservoir(r).MacroNodesID = unique(Reservoir(r).MacroNodesID);
end

%%% ODMACRO
% Macro OD at the city scale. An ODmacro is defined by an origin reservoir
% and a destination reservoir. All existing routes between these reservoirs
% are found in order to assign a particular trip length for all these
% routes in each crossed reservoir.

% ODmacro.OriginID ; ID of the origin reservoir
% ODmacro.DestinationID ; ID of the destination reservoir
% ODmacro.NumPossibleRoutes ; number of all possible routes for the OD
% ODmacro.PossibleRoute ; structure of all possible routes for the OD
% ODmacro.PossibleRoute.ResPath ; list of reservoir IDs, path of the route
% ODmacro.PossibleRoute.NodePath ; list of MacroNode IDs, path of the route
% ODmacro.PossibleRoute.TripLengths ; list of the trip lengths in the reservoirs crossed

% Find origins and destinations (entry or exit through a macro node)
Temp_orires = [];
Temp_orinodes = [];
Temp_destres = [];
Temp_destnodes = [];
for r = 1:NumRes
    for i = Reservoir(r).MacroNodesID
        if strcmp(MacroNode(i).Type,'origin') || strcmp(MacroNode(i).Type,'externalentry') % origins
            Temp_orires = [Temp_orires r];
            Temp_orinodes = [Temp_orinodes i];
        end
        if strcmp(MacroNode(i).Type,'destination') || strcmp(MacroNode(i).Type,'externalexit') % destinations
            Temp_destres = [Temp_destres r];
            Temp_destnodes = [Temp_destnodes i];
        end
    end
end
NumOrigins = length(Temp_orires);
NumDestinations = length(Temp_destres);

NumODmacro = NumOrigins*NumDestinations;
ODmacro = struct('OriginID',cell(1,NumODmacro));
ODmacroID = zeros(NumOrigins,NumDestinations);

od = 1;
for o = 1:NumOrigins % loop on all possible origins
    for d = 1:NumDestinations % loop on all possible destinations
        ODmacro(od).OriginID = o;
        ODmacro(od).DestinationID = d;
        ODmacro(od).ResOriginID = Temp_orires(o);
        ODmacro(od).ResDestinationID = Temp_destres(d);
        ODmacro(od).NodeOriginID = Temp_orinodes(o);
        ODmacro(od).NodeDestinationID = Temp_destnodes(d);
        
        ODmacro(od).NumPossibleRoutes = 0;
        ODmacroID(o,d) = od;
        
        % R1 to R2
        if Temp_orinodes(o) == 1 && Temp_destnodes(d) == 3
            ODmacro(od).NumPossibleRoutes = 5;
            iroute = 1;
            ODmacro(od).PossibleRoute(iroute).ResPath = [1 2 3];
            ODmacro(od).PossibleRoute(iroute).NodePath = [1 5 6 3];
            ODmacro(od).PossibleRoute(iroute).TripLengths = [1000 1000 1000];
            ODmacro(od).PossibleRoute(iroute).NumMicroTrips = 1000;
            iroute = 2;
            ODmacro(od).PossibleRoute(iroute).ResPath = [1 2 3];
            ODmacro(od).PossibleRoute(iroute).NodePath = [1 5 7 3];
            ODmacro(od).PossibleRoute(iroute).TripLengths = [1000 1000 1000];
            ODmacro(od).PossibleRoute(iroute).NumMicroTrips = 1000;
            iroute = 3;
            ODmacro(od).PossibleRoute(iroute).ResPath = [1 2 3];
            ODmacro(od).PossibleRoute(iroute).NodePath = [1 5 6 3];
            ODmacro(od).PossibleRoute(iroute).TripLengths = [1000 1100 1000];
            ODmacro(od).PossibleRoute(iroute).NumMicroTrips = 1000;
            iroute = 4;
            ODmacro(od).PossibleRoute(iroute).ResPath = [1 2 3];
            ODmacro(od).PossibleRoute(iroute).NodePath = [1 5 7 3];
            ODmacro(od).PossibleRoute(iroute).TripLengths = [1000 1200 1000];
            ODmacro(od).PossibleRoute(iroute).NumMicroTrips = 1000;
            iroute = 5;
            ODmacro(od).PossibleRoute(iroute).ResPath = [1 2 3];
            ODmacro(od).PossibleRoute(iroute).NodePath = [1 5 6 3];
            ODmacro(od).PossibleRoute(iroute).TripLengths = [1000 1000 1200];
            ODmacro(od).PossibleRoute(iroute).NumMicroTrips = 1000;
        end
      
        od = od + 1;
    end
end




