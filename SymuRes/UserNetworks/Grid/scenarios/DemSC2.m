%% Scenario 2: loading with a gravity model
%--------------------------------------------------------------------------

Simulation.Duration = 6000; % Simulation duration [s]

Assignment.Periods = [0 Simulation.Duration];
Assignment.PredefRoute = 0;
Assignment.Convergence = 0;

% Distances between reservoirs
Temp_listdist2 = [];
Temp_resdist2 = zeros(NumRes,NumRes);
for i = 1:(NumRes-1)
    for j = (i+1):NumRes
        Temp_d2 = (Reservoir(i).Centroid(1) - Reservoir(j).Centroid(1))^2 + (Reservoir(i).Centroid(2) - Reservoir(j).Centroid(2))^2;
        Temp_listdist2 = [Temp_listdist2 Temp_d2];
        Temp_resdist2(i,j) = Temp_d2;
    end
end
Temp_resdist2 = Temp_resdist2 + Temp_resdist2';
Temp_meandist2 = mean(Temp_listdist2);
Temp_mindist2 = min(Temp_listdist2);
Temp_maxdist2 = max(Temp_listdist2);

% Emission and attraction coefficients
Temp_emicoeffs = zeros(1,NumRes);
Temp_attcoeffs = zeros(1,NumRes);
for i = 1:NumRes
    if (Reservoir(i).RowID == 1 || Reservoir(i).RowID == Nrows) && ...
            (Reservoir(i).ColID == 1 || Reservoir(i).ColID == Ncols)
        % Periphery
        Temp_attcoeffs(i) = 0.2;
        Temp_emicoeffs(i) = 0.8;
    else
        % Center
        Temp_attcoeffs(i) = 0.8;
        Temp_emicoeffs(i) = 0.2;
    end
end

% Peak profile
Td1 = 1000;
Td11 = 2000;
Td12 = 3000;
Td2 = 4000;
q0 = 1;
q1 = 2;
q2 = 1;
qinD = @(t_) q0*trayfct(t_,Td1,Td11,Td12,Td2,q1/q0,q2/q0);

% Demand for all macro OD
%--------------------------------------------------------------------------
for od = 1:NumODmacro
    i = 1;
    ODmacro(od).Demand(i).Purpose = 'cartrip';
    ro = ODmacro(od).ResOriginID;
    rd = ODmacro(od).ResDestinationID;
    Temp_gravcoeff = 4 * Temp_resdist2(ro,rd)/Temp_maxdist2 * (1 - Temp_resdist2(ro,rd)/Temp_maxdist2);
    Temp_dem = 0.04*Temp_emicoeffs(ro)*Temp_attcoeffs(rd)*Temp_gravcoeff;
    ODmacro(od).Demand(i).Time = [0 Td1:60:(Td2+60)]; % [s]
    ODmacro(od).Demand(i).Data = Temp_dem.*[q0 qinD(Td1:60:Td2) q2]; % [veh/s]
end



