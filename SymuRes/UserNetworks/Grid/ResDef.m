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

% Grid model
% Number of grid rows (> 0)
Nrows = 3;
% Number of grid columns (> 0)
Ncols = 3;

% Number of reservoirs
NumRes = Nrows*Ncols;

Reservoir = struct('Centroid',cell(1,NumRes));

% MFD function
% param = [nj nc Pc]
MFDfct = @paraboFD;

% Entry supply function
% param = [nj nc Pc a1*nc a2*nc b*Pc], with 0 < a1 < 1 < a2, and 1 < b
Entryfct = @paraboEntryFD;
% Entryfct = @(n,param) (n <= param(4)).*param(6) + ...
%     (param(4) < n).*(n <= param(5)).*(param(6)+(n-param(4))./(param(5)-param(4)).*(MFDfct(param(5),param(1:3))-param(6))) + ...
%     (param(5) < n).*MFDfct(n,param(1:3));

Exitfct = @paraboExitFD;

% Reservoir definition

% Spacing between reservoirs
ResSpac = 2;

r = 1;
for i = 1:Nrows
    for j = 1:Ncols
        xr = (j-1)*ResSpac;
        yr = 2*ResSpac - (i-1)*ResSpac;
        Reservoir(r).Centroid = [xr yr]; % virtual positioning [x y] for plotting purpose
        Reservoir(r).RowID = i;
        Reservoir(r).ColID = j;
        Reservoir(r).BorderPoints = [xr + ResSpac/2*[-1 1 1 -1 -1]; yr + ResSpac/2*[-1 -1 1 1 -1]];        
        r = r + 1;
    end
end

% All reservoirs
%--------------------------------------------------------------------------
for r = 1:NumRes
    Reservoir(r).AdjacentRes = []; % adjacent reservoirs
    for r2 = 1:NumRes
        Temp_dist = sqrt((Reservoir(r).Centroid(1) - Reservoir(r2).Centroid(1))^2 + (Reservoir(r).Centroid(2) - Reservoir(r2).Centroid(2))^2);
        if Temp_dist < 1.001*ResSpac && r2 ~= r
            Reservoir(r).AdjacentRes = [Reservoir(r).AdjacentRes r2];
        end
    end
    Reservoir(r).FreeflowSpeed = 15; % [m/s]
    Reservoir(r).MaxProd = 3000; % [veh.m/s]
    Reservoir(r).MaxAcc = 1000; % [veh]
    Reservoir(r).CritAcc = 2*Reservoir(r).MaxProd/Reservoir(r).FreeflowSpeed; % for a parabolic MFD
    Reservoir(r).MFDfctParam = [Reservoir(r).MaxAcc Reservoir(r).CritAcc Reservoir(r).MaxProd];
    Reservoir(r).EntryfctParam = [Reservoir(r).MaxAcc Reservoir(r).CritAcc Reservoir(r).MaxProd ...
        0.8*Reservoir(r).CritAcc 1*Reservoir(r).CritAcc 1*Reservoir(r).MaxProd];
    Reservoir(r).ExitfctParam = Reservoir(r).MFDfctParam;
end
if NumRes == 1
    Reservoir(1).AdjacentRes = 1;
end




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

% Origin/entry and destination/exit nodes
i = 1;
for r = 1:NumRes
    MacroNode(i).Type = 'origin';
    MacroNode(i).ResID = r;
    MacroNode(i).Coord = Reservoir(r).Centroid + [0 0.4*ResSpac/2];
    MacroNode(i).Capacity.Time = 0; % [s]
    MacroNode(i).Capacity.Data = 100; % [veh/s]
    i = i + 1;
    
    MacroNode(i).Type = 'destination';
    MacroNode(i).ResID = r;
    MacroNode(i).Coord = Reservoir(r).Centroid + [0 -0.4*ResSpac/2];
    MacroNode(i).Capacity.Time = 0; % [s]
    MacroNode(i).Capacity.Data = 100; % [veh/s]
    i = i + 1;
    
    if Reservoir(r).RowID == 1 || Reservoir(r).RowID == Nrows
        if Reservoir(r).RowID == 1
            Temp_dir = 1;
        elseif Reservoir(r).RowID == Nrows
            Temp_dir = -1;
        end
        MacroNode(i).Type = 'externalentry';
        MacroNode(i).ResID = r;
        MacroNode(i).Coord = Reservoir(r).Centroid + [-0.2*ResSpac/2 Temp_dir*ResSpac/2];
        MacroNode(i).Capacity.Time = 0; % [s]
        MacroNode(i).Capacity.Data = 100; % [veh/s]
        i = i + 1;
        MacroNode(i).Type = 'externalexit';
        MacroNode(i).ResID = r;
        MacroNode(i).Coord = Reservoir(r).Centroid + [0.2*ResSpac/2 Temp_dir*ResSpac/2];
        MacroNode(i).Capacity.Time = 0; % [s]
        MacroNode(i).Capacity.Data = 100; % [veh/s]
        i = i + 1;
    end
    
    if Reservoir(r).ColID == 1 || Reservoir(r).ColID == Ncols
        if Reservoir(r).ColID == 1
            Temp_dir = -1;
        elseif Reservoir(r).ColID == Ncols
            Temp_dir = 1;
        end
        MacroNode(i).Type = 'externalentry';
        MacroNode(i).ResID = r;
        MacroNode(i).Coord = Reservoir(r).Centroid + [Temp_dir*ResSpac/2 0.2*ResSpac/2];
        MacroNode(i).Capacity.Time = 0; % [s]
        MacroNode(i).Capacity.Data = 100; % [veh/s]
        i = i + 1;
        MacroNode(i).Type = 'externalexit';
        MacroNode(i).ResID = r;
        MacroNode(i).Coord = Reservoir(r).Centroid + [Temp_dir*ResSpac/2 -0.2*ResSpac/2];
        MacroNode(i).Capacity.Time = 0; % [s]
        MacroNode(i).Capacity.Data = 100; % [veh/s]
        i = i + 1;
    end
end

% Border nodes
for r = 1:(NumRes-1)
    for r2 = (r+1):NumRes
        if ~isempty(find(Reservoir(r).AdjacentRes == r2,1))
            xr = Reservoir(r).Centroid(1);
            yr = Reservoir(r).Centroid(2);
            xr2 = Reservoir(r2).Centroid(1);
            yr2 = Reservoir(r2).Centroid(2);
            
            % From r to r2
            MacroNode(i).Type = 'border';
            MacroNode(i).ResID = [r r2];
            MacroNode(i).Coord = [(xr + xr2)/2 + 0.2*(yr - yr2) (yr + yr2)/2 + 0.2*(xr2 - xr)];
            MacroNode(i).Capacity.Time = 0; % [s]
            MacroNode(i).Capacity.Data = 100; % [veh/s]
            i = i + 1;
            
            % From r2 to r
            MacroNode(i).Type = 'border';
            MacroNode(i).ResID = [r2 r];
            MacroNode(i).Coord = [(xr + xr2)/2 - 0.2*(yr - yr2) (yr + yr2)/2 - 0.2*(xr2 - xr)];
            MacroNode(i).Capacity.Time = 0; % [s]
            MacroNode(i).Capacity.Data = 100; % [veh/s]
            i = i + 1;
        end
    end
end

NumMacroNodes = i - 1;

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

Temp_CostMatrix = Inf*ones(NumMacroNodes,NumMacroNodes);
for i = 1:NumMacroNodes
    if strcmp(MacroNode(i).Type,'border')
        for j = 1:NumMacroNodes
            if MacroNode(i).ResID(2) == MacroNode(j).ResID(1)
                if ~strcmp(MacroNode(j).Type,'origin') && ~strcmp(MacroNode(j).Type,'externalentry')
                    Temp_CostMatrix(i,j) = 1;
                end
            end
        end
    elseif strcmp(MacroNode(i).Type,'origin') || strcmp(MacroNode(i).Type,'externalentry')
        for j = 1:NumMacroNodes
            if MacroNode(i).ResID == MacroNode(j).ResID(1)
                if ~strcmp(MacroNode(j).Type,'origin') && ~strcmp(MacroNode(j).Type,'externalentry')
                    Temp_CostMatrix(i,j) = 1;
                end
            end
        end
    end
end

od = 1;
for o = 1:NumOrigins % loop on all possible origins
    for d = 1:NumDestinations % loop on all possible destinations
        ODmacro(od).OriginID = o;
        ODmacro(od).DestinationID = d;
        ODmacro(od).ResOriginID = Temp_orires(o);
        ODmacro(od).ResDestinationID = Temp_destres(d);
        ODmacro(od).NodeOriginID = Temp_orinodes(o);
        ODmacro(od).NodeDestinationID = Temp_destnodes(d);
        ODmacroID(o,d) = od;
        
        Temp_shortestPaths = kShortestPath(Temp_CostMatrix,Temp_orinodes(o),Temp_destnodes(d),3);
        iroute = 1;
        for ipath = 1:length(Temp_shortestPaths)
            Temp_nodepath = Temp_shortestPaths{ipath};
            if length(Temp_nodepath) <= 2
                Temp_respath = MacroNode(Temp_nodepath(1)).ResID(1);
                ODmacro(od).PossibleRoute(iroute).ResPath = Temp_respath;
                ODmacro(od).PossibleRoute(iroute).NodePath = Temp_nodepath;
                ODmacro(od).PossibleRoute(iroute).TripLengths = 1000*ones(1,length(Temp_respath));
                ODmacro(od).PossibleRoute(iroute).NumMicroTrips = 1000;
                iroute = iroute + 1;
            else
                Temp_respath = [];
                for i = Temp_nodepath(2:end)
                    Temp_respath = [Temp_respath MacroNode(i).ResID(1)];
                end
                if isequal(Temp_respath,unique(Temp_respath,'stable'))
                    ODmacro(od).PossibleRoute(iroute).ResPath = Temp_respath;
                    ODmacro(od).PossibleRoute(iroute).NodePath = Temp_nodepath;
                    ODmacro(od).PossibleRoute(iroute).TripLengths = 1000*ones(1,length(Temp_respath));
                    ODmacro(od).PossibleRoute(iroute).NumMicroTrips = 1000;
                    iroute = iroute + 1;
                end
            end
        end
        ODmacro(od).NumPossibleRoutes = iroute - 1;
        
        od = od + 1;
    end
end




