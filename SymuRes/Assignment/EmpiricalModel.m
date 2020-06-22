function ODflowdistrib = EmpiricalModel(od,RoutesList,Reservoir,Route,Assignment,Simulation)
% ODflowdistrib = EmpiricalModel(od,RoutesList,Reservoir,Route,Assignment,Simulation)
% Return the path flow distribution for a given OD according to a traffic
% equilibrium model
%
% Nov 2019 - Guilhem Mariotte
%
% INPUTS
%---- od         : scalar, ID of the considered OD
%---- RoutesList : vector, list of route IDs for the considered OD
%---- Reservoir  : Reservoir structure
%---- Route      : Route structure
%---- Assignment : Assignment structure
%---- Simulation : Simulation structure
%
% OUTPUTS
%---- ODflowdistrib : vector, OD path flow distribution, i.e. assignment
%                     coefficient of the routes considered

Temp_RouteIDs = RoutesList;
Temp_Nroutes = length(RoutesList);

if Assignment.model == 100
    % Manual assignment
    ODflowdistrib = Assignment.ManualCoefficients;
    
elseif Assignment.model == 101
    % Equi-probability
    ODflowdistrib(Temp_RouteIDs) = ones(1,Temp_Nroutes)./Temp_Nroutes;
    
elseif Assignment.model == 102
    % Prorata of the number of micro trips
    Temp_Ntrips = zeros(1,Temp_Nroutes);
    i_r = 1;
    for iroute = Temp_RouteIDs
        Temp_Ntrips(i_r) = Route(iroute).NumMicroTrips;
        i_r = i_r + 1;
    end
    ODflowdistrib(Temp_RouteIDs) = Temp_Ntrips./sum(Temp_Ntrips);
    
elseif Assignment.model == 103
    % Weighting on reservoirs
    Assignment.model = 1; % Apply Wardrop's principle
    Assignment.Behavior = 1;
    ODflowdistrib = ChoiceModel(od,Temp_RouteIDs,Reservoir,Route,Assignment,Simulation);
    Temp_weight = ODflowdistrib;
    Temp_totweight = 1;
    for iroute = Temp_RouteIDs
        for r = Route(iroute).ResPath
            Temp_weight(iroute) = Temp_weight(iroute) + Assignment.ReservoirWeight(r);
            Temp_totweight = Temp_totweight + Assignment.ReservoirWeight(r);
        end
    end
    ODflowdistrib = Temp_weight./Temp_totweight; % bonus when crossing reservoir, according to their weight parameter
    
elseif Assignment.model == 104
    % Penalty on a routes
    Assignment.model = 1; % Apply Wardrop's principle
    Assignment.Behavior = 1;
    ODflowdistrib = ChoiceModel(od,Temp_RouteIDs,Reservoir,Route,Assignment,Simulation);
    if length(Temp_RouteIDs) > 1
        Temp_thres = Assignment.PenalizedAssignCoeff;
        Temp_coeff = ODflowdistrib(Temp_RouteIDs);
        Temp_coefftot = 1;
        Temp_indexes = [];
        i_r = 1;
        for iroute = Temp_RouteIDs
            Temp_haspenalty = 0;
            for i_p = 1:length(Assignment.PenalizedResPath)
                if isequal(Route(iroute).ResPath,Assignment.PenalizedResPath{i_p}) && ODflowdistrib(iroute) > Temp_thres
                    Temp_haspenalty = 1;
                    Temp_coeff(i_r) = Temp_thres; % apply a penalty on some routes
                    Temp_coefftot = Temp_coefftot - Temp_coeff(i_r);
                end
            end
            if Temp_haspenalty == 0
                Temp_indexes = [Temp_indexes i_r];
            end
            i_r = i_r + 1;
        end
        if length(Temp_indexes) < length(Temp_RouteIDs)
            Temp_flowdistrib = ChoiceModel(od,Temp_RouteIDs(Temp_indexes),Reservoir,Route,Assignment,Simulation);
            Temp_coeff(Temp_indexes) = Temp_flowdistrib(Temp_RouteIDs(Temp_indexes));
            Temp_coeff(Temp_indexes) = Temp_coeff(Temp_indexes).*Temp_coefftot;
            ODflowdistrib(Temp_RouteIDs) = Temp_coeff;
        end
    end
    
end

end
