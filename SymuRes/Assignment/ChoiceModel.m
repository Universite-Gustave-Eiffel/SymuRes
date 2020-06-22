function ODflowdistrib = ChoiceModel(od,RoutesList,Reservoir,Route,Assignment,Simulation)
% ODflowdistrib = ChoiceModel(od,RoutesList,Reservoir,Route,Assignment,Simulation)
% Return the path flow distribution for a given OD according to a traffic
% equilibrium model
%
% Feb 2019 - Sergio F. A. Batista
%
% Reference:
% Batista & Leclercq (TS, 2019)
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

Temp_iter = Assignment.CurIteration;
Temp_RouteIDs = RoutesList;

% Assignment period time window
% if Temp_iter > 1
%     NumTimes = length(Simulation.CurrentTime);
%     Temp_StartTimeID = floor(Assignment.CurrentTime/Simulation.TimeStep) + 1;
%     Temp_EndTimeID = min([floor((Assignment.CurrentTime + Assignment.Period)/Simulation.TimeStep) NumTimes-1]);
% end

% Assignment period time window
NumTimes = length(Simulation.Time);
Temp_StartTimeID = floor(Assignment.CurrentTime/Simulation.TimeStep) + 1;
Temp_EndTimeID = min([floor(Assignment.Periods(Assignment.CurrentPeriodID+1)/Simulation.TimeStep) NumTimes-1]);
% Previous assignment period time window
if Temp_StartTimeID == 1
    Temp_StartTimeID_Previous = Temp_StartTimeID;
    Temp_EndTimeID_Previous = Temp_EndTimeID;
elseif Temp_StartTimeID > 1
    Temp_StartTimeID_Previous = Temp_StartTimeID - (Temp_EndTimeID - Temp_StartTimeID) - 1;
    Temp_EndTimeID_Previous = Temp_StartTimeID - 1;
end

DistributionPerRoute = zeros(1,length(Route));


%% Models from the literature
%--------------------------------------------------------------------------

% Update the AL for the bounded rationality:
if Assignment.Behavior==2
    Utility=zeros(1,length(Temp_RouteIDs));
    for k=1:length(Temp_RouteIDs)
        macro_path=Route(Temp_RouteIDs(k)).ResPath;
        path_length=Route(Temp_RouteIDs(k)).TripLengths;
        path_speed=[];
        for i_res=1:length(macro_path)
            if Temp_iter == 1
                if Temp_StartTimeID == 1
                    Temp_meanV = Reservoir(macro_path(i_res)).FreeflowSpeed;
                else
                    Temp_meanV = mean(Reservoir(macro_path(i_res)).MeanSpeed(Temp_StartTimeID_Previous:Temp_EndTimeID_Previous));
                end
            elseif Temp_iter > 1
                Temp_meanV = mean(Reservoir(macro_path(i_res)).MeanSpeed(Temp_StartTimeID:Temp_EndTimeID));
            end
            if Temp_meanV > 0
                path_speed=[path_speed Temp_meanV];
            else
                path_speed=[path_speed Reservoir(macro_path(i_res)).FreeflowSpeed];
            end
        end
        Utility(k)=sum(path_length./path_speed);
    end
    
    % Updating the AL for the od pair.
    Assignment.AL(od)=min(Utility)+Assignment.ALparam.*min(Utility);
    
    clear Utility macro_path k path_length path_speed i_res
end


% Assignment models:
% Deterministic User Equilibrium (DUE) based on the first Wardrop
% principle. Or, the BR-DUE based on a random search order among
% the satisficing routes.
%----------------------------------------------------------------
if Assignment.model==1
    
    % Calculating the utility function.
    Utility=zeros(1,length(Temp_RouteIDs));
    
    for k=1:length(Temp_RouteIDs)
        macro_path=Route(Temp_RouteIDs(k)).ResPath;
        path_length=Route(Temp_RouteIDs(k)).TripLengths;
        path_speed=[];
        for i_res=1:length(macro_path)
            if Temp_iter == 1
                if Temp_StartTimeID == 1
                    Temp_meanV = Reservoir(macro_path(i_res)).FreeflowSpeed;
                else
                    Temp_meanV = mean(Reservoir(macro_path(i_res)).MeanSpeed(Temp_StartTimeID_Previous:Temp_EndTimeID_Previous));
                end
            elseif Temp_iter > 1
                Temp_meanV = mean(Reservoir(macro_path(i_res)).MeanSpeed(Temp_StartTimeID:Temp_EndTimeID));
            end
            if Temp_meanV > 0
                path_speed=[path_speed Temp_meanV];
            else
                path_speed=[path_speed Reservoir(macro_path(i_res)).FreeflowSpeed];
            end
        end
        Utility(k)=sum(path_length./path_speed);
    end
    
    % DUE:
    if Assignment.Behavior==1
        for k=1:length(Temp_RouteIDs)
            if (Utility(k)-min(Utility))<0.1 % To avoid numerical fluctuations
                DistributionPerRoute(Temp_RouteIDs(k)) = 1/length(find(Utility - min(Utility) < 0.1));
            else
                DistributionPerRoute(Temp_RouteIDs(k)) = 0;
            end
        end
        
        % BR-DUE
    elseif Assignment.Behavior==2
        
        Satisficing=[];
        % Update the utilities and calculating the satisficing alternatives:
        for k=1:length(Temp_RouteIDs)
            if Utility(k)<=Assignment.AL(od)+0.01
                Satisficing=[Satisficing Temp_RouteIDs(k)];
            end
        end
        
        Choices=zeros(Assignment.parameters(4),1);
        for i_choices=1:Assignment.parameters(4)
            
            % If there are satisficing routes
            if sum(Satisficing)>0
                
                % Constructing the ranges:
                range=[];
                range=[range, 0];
                for k=1:length(Satisficing)
                    range = [range, k/length(Satisficing)];
                end
                
                % Choosing the alternative
                r = rand();
                m1 = (r >= range);
                m2 = (r <= range);
                m2(1:end-1) = m2(2:end);
                m2(end) = 0;
                
                Choices(i_choices)=Satisficing(find(m1 & m2));
                
                clear range m1 m2
                
                % If there are no satisficing alternatives, the users choose the best alternative available. That is,
                % the minimum utility route.
            else
                [aa,bb] = min(Utility(i_samples,:));
                Choices(i_choices)=Temp_RouteIDs(bb);
            end
            
        end
        
        for i_routes=1:length(Temp_RouteIDs)
            DistributionPerRoute(Temp_RouteIDs(i_routes))=sum(Choices(:)==i_routes)/Assignment.parameters(4);
        end
        
        % RT-DUE
    elseif Assignment.Behavior==3
        
        % Perceived regret.
        Regret=zeros(1,length(Temp_RouteIDs));
        
        % Calculating the perceived regret:
        Regret=Utility+exp(-Assignment.RTparam.*Utility).*(Utility-min(Utility));
        
        % Assign users based on the minimization regret behavior:
        for k=1:length(Temp_RouteIDs)
            if (Regret(k)-min(Regret))<0.1 % To avoid numerical fluctuations
                DistributionPerRoute(Temp_RouteIDs(k)) = 1/length(find(Regret - min(Regret) < 0.1));
            else
                DistributionPerRoute(Temp_RouteIDs(k)) = 0;
            end
        end
        
    end
    
    clear Utility k macro_path path_length path_speed i_res Satisficing range r m1 m2 Choices aa bb i_choices Regret
end

% Stochastic User Equilibrium that is solved using Monte Carlo
% simulations.
%-------------------------------------------------------------
if Assignment.model==2
    
    % The first speed consists of gathering the data samples,
    % depending on the utility function chosen by the user.
    % Assignment.parameters(1)=1 - Uncertainty on the trip lengths.
    % Assignment.parameters(1)=2 - Uncertainty on the mean speed.
    % Assignment.parameters(1)=3 - Uncertainty on the trip lengths and mean speed.
    
    Utility=zeros(Assignment.parameters(3),length(Temp_RouteIDs));
    
    if Assignment.Behavior==1 || Assignment.Behavior==3
        Route_choice=zeros(Assignment.parameters(2),1);
    elseif Assignment.Behavior==2
        Route_choice=zeros(Assignment.parameters(4),length(Temp_RouteIDs));
    end
    
    for k=1:length(Temp_RouteIDs)
        
        macro_path=Route(Temp_RouteIDs(k)).ResPath;
        
        if Assignment.UtilityID==1 || Assignment.UtilityID==3
            TripLength_samples=zeros(Assignment.parameters(1),length(macro_path));
        end
        if Assignment.UtilityID==2 || Assignment.UtilityID==3
            MeanSpeed_samples=zeros(Assignment.parameters(2),length(macro_path));
        end
        
        % Gathering the Mean Speed inside the reservoirs.
        path_speed=[];
        for i_res=1:length(macro_path)
            if Temp_iter == 1
                if Temp_StartTimeID == 1
                    Temp_meanV = Reservoir(macro_path(i_res)).FreeflowSpeed;
                else
                    Temp_meanV = mean(Reservoir(macro_path(i_res)).MeanSpeed(Temp_StartTimeID_Previous:Temp_EndTimeID_Previous));
                end
            elseif Temp_iter > 1
                Temp_meanV = mean(Reservoir(macro_path(i_res)).MeanSpeed(Temp_StartTimeID:Temp_EndTimeID));
            end
            if Temp_meanV > 0
                path_speed=[path_speed Temp_meanV];
            else
                path_speed=[path_speed Reservoir(macro_path(i_res)).FreeflowSpeed];
            end
        end
        
        % Gathering the trip lengths inside the reservoirs:
        path_length=Route(Temp_RouteIDs(k)).TripLengths;
        
        % Gathering the samples of the mean speed inside the reservoirs:
        if Assignment.UtilityID==2 || Assignment.UtilityID==3
            for i_res=1:length(macro_path)
                if Temp_StartTimeID == 1
                    MeanSpeed_samples(:,i_res)=Reservoir(macro_path(i_res)).FreeflowSpeed;
                elseif Temp_StartTimeID > 1
                    if Temp_iter==1
                        for i_samples=1:Assignment.parameters(2)
                            r=Temp_StartTimeID_Previous+round(rand()*length(Reservoir(macro_path(i_res)).MeanSpeed(Temp_StartTimeID_Previous:Temp_EndTimeID_Previous))+0.5);
                            if r>Temp_EndTimeID_Previous
                                r=Temp_EndTimeID_Previous;
                            elseif r<Temp_StartTimeID_Previous
                                r=Temp_StartTimeID_Previous;
                            end
                            MeanSpeed_samples(i_samples,i_res)=Reservoir(macro_path(i_res)).MeanSpeed(r);
                        end
                    elseif Temp_iter>1
                        for i_samples=1:Assignment.parameters(2)
                            r=Temp_StartTimeID+round(rand()*length(Reservoir(macro_path(i_res)).MeanSpeed(Temp_StartTimeID:Temp_EndTimeID))+0.5);
                            if r>Temp_EndTimeID
                                r=Temp_EndTimeID;
                            elseif r<Temp_StartTimeID
                                r=Temp_StartTimeID;
                            end
                            MeanSpeed_samples(i_samples,i_res)=Reservoir(macro_path(i_res)).MeanSpeed(r);
                        end
                    end
                end
            end
        end
        
        % Getting the samples of the trip lengths:
        if Assignment.UtilityID==1 || Assignment.UtilityID==3
            if Assignment.TripLengthMethod==1 % Method 1.
                for i_res=1:length(macro_path)
                    for i_samples=1:Assignment.parameters(1)
                        r=round(rand()*length(TripLengths_per_reservoir_1(macro_path(i_res)).Length)+0.5);
                        TripLength_samples(i_samples,i_res)=TripLengths_per_reservoir_1(macro_path(i_res)).Length{r};
                    end
                end
            elseif Assignment.TripLengthMethod==2 % Method 2.
                if length(macro_path)==1
                    for i_loop=1:length(TripLengths_per_reservoir_per_bin_2)
                        if macro_path(1)==TripLengths_per_reservoir_per_bin_2(i_loop).CurrentRes && macro_path(1)==TripLengths_per_reservoir_per_bin_2(i_loop).DestinationRes
                            for i_samples=1:Assignment.parameters(1)
                                r=round(rand()*length(TripLengths_per_reservoir_per_bin_2(i_loop).Length)+0.5);
                                TripLength_samples(i_samples,1)=TripLengths_per_reservoir_per_bin_2(i_loop).Length{r};
                            end
                        end
                    end
                elseif length(macro_path)>1
                    % Intermediate reservoirs:
                    for i_res=1:length(macro_path)-1
                        for i_loop=1:length(TripLengths_per_reservoir_per_bin_2)
                            if macro_path(i_res)==TripLengths_per_reservoir_per_bin_2(i_loop).CurrentRes && macro_path(i_res+1)==TripLengths_per_reservoir_per_bin_2(i_loop).DestinationRes
                                for i_samples=1:Assignment.parameters(1)
                                    r=round(rand()*length(TripLengths_per_reservoir_per_bin_2(i_loop).Length)+0.5);
                                    TripLength_samples(i_samples,i_res)=TripLengths_per_reservoir_per_bin_2(i_loop).Length{r};
                                end
                            end
                        end
                    end
                    % Destination reservoirs of the macro-path:
                    for i_loop=1:length(TripLengths_per_reservoir_per_bin_2)
                        if macro_path(length(macro_path))==TripLengths_per_reservoir_per_bin_2(i_loop).CurrentRes && macro_path(length(macro_path))==TripLengths_per_reservoir_per_bin_2(i_loop).DestinationRes
                            for i_samples=1:Assignment.parameters(1)
                                r=round(rand()*length(TripLengths_per_reservoir_per_bin_2(i_loop).Length)+0.5);
                                TripLength_samples(i_samples,i_res+1)=TripLengths_per_reservoir_per_bin_2(i_loop).Length{r};
                            end
                        end
                    end
                end
            elseif Assignment.TripLengthMethod==3 % Method 3.
                if length(macro_path)==1
                    for i_loop=counter1+counter2+counter3+1:length(TripLengths_per_origin_per_bin_per_destination_3)
                        if TripLengths_per_origin_per_bin_per_destination_3(i_loop).i_res==macro_path(1)
                            for i_samples=1:Assignment.parameters(1)
                                r=round(rand()*length(TripLengths_per_origin_per_bin_per_destination_3(i_loop).Length)+0.5);
                                TripLength_samples(i_samples,1)=TripLengths_per_origin_per_bin_per_destination_3(i_loop).Length(r);
                            end
                        end
                    end
                elseif length(macro_path)==2
                    % Origin reservoir to the border:
                    for i_loop=counter1+counter2+1:counter1+counter2+counter3
                        if TripLengths_per_origin_per_bin_per_destination_3(i_loop).h_res==macro_path(1) && TripLengths_per_origin_per_bin_per_destination_3(i_loop).i_res==macro_path(2)
                            for i_samples=1:Assignment.parameters(1)
                                r=round(rand()*length(TripLengths_per_origin_per_bin_per_destination_3(i_loop).Length)+0.5);
                                TripLength_samples(i_samples,1)=TripLengths_per_origin_per_bin_per_destination_3(i_loop).Length(r);
                            end
                        end
                    end
                    % Border to destination reservoir.
                    for i_loop=counter1+1:counter1+counter2
                        if TripLengths_per_origin_per_bin_per_destination_3(i_loop).i_res==macro_path(2) && TripLengths_per_origin_per_bin_per_destination_3(i_loop).h_res==macro_path(1)
                            for i_samples=1:Assignment.parameters(1)
                                r=round(rand()*length(TripLengths_per_origin_per_bin_per_destination_3(i_loop).Length)+0.5);
                                TripLength_samples(i_samples,2)=TripLengths_per_origin_per_bin_per_destination_3(i_loop).Length(r);
                            end
                        end
                    end
                elseif length(macro_path)>2
                    % Origin reservoir to the border:
                    for i_loop=counter1+counter2+1:counter1+counter2+counter3
                        if TripLengths_per_origin_per_bin_per_destination_3(i_loop).h_res==macro_path(1) && TripLengths_per_origin_per_bin_per_destination_3(i_loop).i_res==macro_path(2)
                            for i_samples=1:Assignment.parameters(1)
                                r=round(rand()*length(TripLengths_per_origin_per_bin_per_destination_3(i_loop).Length)+0.5);
                                TripLength_samples(i_samples,1)=TripLengths_per_origin_per_bin_per_destination_3(i_loop).Length(r);
                            end
                        end
                    end
                    for i_res = 1:length(macro_path)-2
                        for i_loop=1:counter1
                            if TripLengths_per_origin_per_bin_per_destination_3(i_loop).h_res==macro_path(i_res) && TripLengths_per_origin_per_bin_per_destination_3(i_loop).i_res==macro_path(i_res+1) && TripLengths_per_origin_per_bin_per_destination_3(i_loop).j_res==macro_path(i_res+2)
                                for i_samples=1:Assignment.parameters(1)
                                    r=round(rand()*length(TripLengths_per_origin_per_bin_per_destination_3(i_loop).Length)+0.5);
                                    TripLength_samples(i_samples,i_res+1)=TripLengths_per_origin_per_bin_per_destination_3(i_loop).Length(r);
                                end
                            end
                        end
                    end
                    % Border to destination reservoir.
                    for i_loop=counter1+1:counter1+counter2
                        if TripLengths_per_origin_per_bin_per_destination_3(i_loop).i_res==macro_path(length(macro_path)) && TripLengths_per_origin_per_bin_per_destination_3(i_loop).h_res==macro_path(length(macro_path)-1)
                            for i_samples=1:Assignment.parameters(1)
                                r=round(rand()*length(TripLengths_per_origin_per_bin_per_destination_3(i_loop).Length)+0.5);
                                TripLength_samples(i_samples,length(macro_path))=TripLengths_per_origin_per_bin_per_destination_3(i_loop).Length(r);
                            end
                        end
                    end
                end
            elseif Assignment.TripLengthMethod==4 % Method 4.
                for i_path=1:length(TripLengths_per_MacroPath_4)
                    if isequal(macro_path,TripLengths_per_MacroPath_4(i_path).MacroPath)==1
                        for i_res=1:length(macro_path)
                            for i_samples=1:Assignment.parameters(1)
                                r=round(rand()*TripLengths_per_MacroPath_4(i_path).Ntrips);
                                if r==0
                                    r=1;
                                end
                                TripLength_samples(i_samples,i_res)=TripLengths_per_MacroPath_4(i_path).Length{i_res}(r);
                            end
                        end
                    end
                end
            end
        end
        
        if Assignment.UtilityID==1
            for i_samples=1:Assignment.parameters(3)
                r1=round(rand()*Assignment.parameters(1)+0.5);
                Utility(i_samples,k)=sum(TripLength_samples(r1,:)./path_speed);
            end
        elseif Assignment.UtilityID==2
            for i_samples=1:Assignment.parameters(3)
                r2=round(rand()*Assignment.parameters(2)+0.5);
                Utility(i_samples,k)=sum((path_length.*MeanSpeed_samples(r2,:))./(path_speed.^2));
            end
        elseif Assignment.UtilityID==3
            for i_samples=1:Assignment.parameters(3)
                r1=round(rand()*Assignment.parameters(1)+0.5);
                r2=round(rand()*Assignment.parameters(2)+0.5);
                Utility(i_samples,k)=sum(TripLength_samples(r1,:)./path_speed) + sum((path_length.*MeanSpeed_samples(r2,:))./(path_speed.^2));
            end
        end
    end
    
    % SUE:
    if Assignment.Behavior==1
        for i_samples=1:length(Utility)
            [aa,Route_choice(i_samples)] = min(Utility(i_samples,:));
        end
        
        % BR-SUE
    elseif Assignment.Behavior==2
        
        for i_samples=1:length(Utility)
            
            Satisficing=[];
            % Update the utilities and calculating the satisficing alternatives:
            for k=1:length(Temp_RouteIDs)
                if Utility(k)<=Assignment.AL(od)+0.01
                    Satisficing=[Satisficing Temp_RouteIDs(k)];
                end
            end
            
            Choices=zeros(Assignment.parameters(4),1);
            for i_choices=1:Assignment.parameters(4)
                
                % If there are satisficing routes
                if sum(Satisficing)>0
                    
                    % Constructing the ranges:
                    range=[];
                    range=[range, 0];
                    for k=1:length(Satisficing)
                        range = [range, k/length(Satisficing)];
                    end
                    
                    % Choosing the alternative
                    r = rand();
                    m1 = (r >= range);
                    m2 = (r <= range);
                    m2(1:end-1) = m2(2:end);
                    m2(end) = 0;
                    
                    Choices(i_choices)=Satisficing(find(m1 & m2));
                    
                    clear range m1 m2
                    
                    % If there are no satisficing alternatives, the users choose the best alternative available. That is,
                    % the minimum utility route.
                else
                    [aa,bb] = min(Utility(i_samples,:));
                    Choices(i_choices)=Temp_RouteIDs(bb);
                end
            end
            
            for i_route=1:length(Temp_RouteIDs)
                Route_choice(i_samples,i_route)=sum(Choices(:)==Temp_RouteIDs(i_route))/Assignment.parameters(4);
            end
            
        end
        
        % RT-SUE
    elseif Assignment.Behavior==3
        
        % Perceived regret.
        Regret=zeros(Assignment.parameters(3),length(Temp_RouteIDs));
        
        % Calculating the utility function.
        TT_mean=zeros(1,length(Temp_RouteIDs));
        
        for k=1:length(Temp_RouteIDs)
            macro_path=Route(Temp_RouteIDs(k)).ResPath;
            path_length=Route(Temp_RouteIDs(k)).TripLengths;
            path_speed=[];
            for i_res=1:length(macro_path)
                if Temp_iter == 1
                    if Temp_StartTimeID == 1
                        Temp_meanV = Reservoir(macro_path(i_res)).FreeflowSpeed;
                    else
                        Temp_meanV = mean(Reservoir(macro_path(i_res)).MeanSpeed(Temp_StartTimeID_Previous:Temp_EndTimeID_Previous));
                    end
                elseif Temp_iter > 1
                    Temp_meanV = mean(Reservoir(macro_path(i_res)).MeanSpeed(Temp_StartTimeID:Temp_EndTimeID));
                end
                if Temp_meanV > 0
                    path_speed=[path_speed Temp_meanV];
                else
                    path_speed=[path_speed Reservoir(macro_path(i_res)).FreeflowSpeed];
                end
            end
            TT_mean(k)=sum(path_length./path_speed);
        end
        clear k macro_path path_length path_speed i_res Temp_meanV
        
        aux=min(TT_mean)-TT_mean;
        for k=1:length(aux)
            if aux(k)==0
                aux(k)=1E-15;
            end
        end
        
        % Calculating the perceived regret:
        for i_samples=1:Assignment.parameters(3)
            Regret(i_samples,:)=Utility(i_samples,:)-(1-exp(-Assignment.RTparam.*aux));
        end
        
        % Assign users based on the minimization regret behavior:
        for i_samples=1:length(Regret)
            [aa,Route_choice(i_samples)] = min(Regret(i_samples,:));
        end
        
    end
    clear i_choices Choices i_route Choices aa bb i_choices i_samples Regret
    
    % Calculate the new assignment flows:
    if Assignment.Behavior==1 || Assignment.Behavior==3
        for i_routes=1:length(Temp_RouteIDs)
            DistributionPerRoute(Temp_RouteIDs(i_routes))=sum(Route_choice(:)==i_routes)/Assignment.parameters(3);
        end
    elseif Assignment.Behavior==2
        for i_routes=1:length(Temp_RouteIDs)
            DistributionPerRoute(Temp_RouteIDs(i_routes))=sum(Route_choice(:,i_routes))/Assignment.parameters(3);
        end
    end
    
    clear Utility Route_choice k macro_path TripLength_samples path_speed i_res i_samples r i_loop i_res i_path MeanSpeed_samples path_length
    clear r1 r2 aa Satisficing i_routes
    
end

% Mean-variance model:
%---------------------
if Assignment.model==3
    
    Utility=zeros(1,length(Temp_RouteIDs));
    for k=1:length(Temp_RouteIDs)
        macro_path=Route(Temp_RouteIDs(k)).ResPath;
        
        % Recover the lengths:
        for i_path=1:length(TripLengths_per_MacroPath_4)
            if isequal(macro_path,TripLengths_per_MacroPath_4(i_path).MacroPath)==1
                auxlengths = zeros(length(TripLengths_per_MacroPath_4(i_path).Length),length(macro_path));
                for j_res=1:length(macro_path)
                    auxlengths(:,j_res)=TripLengths_per_MacroPath_4(i_path).Length(j_res,:);
                end
            end
        end
        
        % Recovering the speeds:
        auxspeeds = zeros(Temp_EndTimeID-Temp_StartTimeID+1,length(macro_path));
        for i_res=1:length(macro_path)
            if Temp_StartTimeID == 1
                auxspeeds(:,i_res)=Reservoir(macro_path(i_res)).FreeflowSpeed;
            else
                if Temp_iter==1
                    auxspeeds(:,i_res)=Reservoir(macro_path(i_res)).MeanSpeed(Temp_StartTimeID_Previous:Temp_EndTimeID_Previous);
                elseif Temp_iter>1
                    auxspeeds(:,i_res)=Reservoir(macro_path(i_res)).MeanSpeed(Temp_StartTimeID:Temp_EndTimeID);
                end
            end
        end
        
        covariance_vector = zeros(1,length(macro_path));
        pathlegnth_mean = zeros(1,length(macro_path));
        pathspeed_mean = zeros(1,length(macro_path));
        pathlegnth_var = zeros(1,length(macro_path));
        pathspeed_var = zeros(1,length(macro_path));
        for i_res=1:length(macro_path)
            pathspeeds = zeros(1,Assignment.parameters(1));
            pathlength = zeros(1,Assignment.parameters(1));
            for i_samples=1:Assignment.parameters(1)
                pathlength(1,i_samples)=auxlengths(round(rand()*length(auxlengths)+0.5),i_res);
                pathspeeds(1,i_samples)=auxspeeds(round(rand()*length(auxspeeds)+0.5),i_res);
            end
            a = cov(pathlength(1,:),pathspeeds(1,:));
            pathlegnth_var(1,i_res) = a(1,1);
            pathspeed_var(1,i_res) = a(2,2);
            covariance_vector(1,i_res) = a(1,2);
            clear a
            pathlegnth_mean(1,i_res) = mean(pathlength(1,:));
            pathspeed_mean(1,i_res) = mean(pathspeeds(1,:));
        end
        clear pathlength pathspeeds
        
        Utility(k)=Assignment.VOT*sum((pathlegnth_mean./pathspeed_mean))+Assignment.VOR*(sum((pathlegnth_mean./pathspeed_mean).^2)*(sum(pathlegnth_var./pathlegnth_mean.^2)+sum(pathspeed_var./pathspeed_mean.^2)+sum(2.*covariance_vector./(pathlegnth_mean.*pathspeed_mean))));
    end
    
    % DUE:
    if Assignment.Behavior==1
        for k=1:length(Temp_RouteIDs)
            if (Utility(k)-min(Utility))<0.1 % To avoid numerical fluctuations
                DistributionPerRoute(Temp_RouteIDs(k)) = 1/length(find(Utility - min(Utility) < 0.1));
            else
                DistributionPerRoute(Temp_RouteIDs(k)) = 0;
            end
        end
        
        % BR-DUE
    elseif Assignment.Behavior==2
        
        Satisficing=[];
        % Update the utilities and calculating the satisficing alternatives:
        for k=1:length(Temp_RouteIDs)
            if Utility(k)<=Assignment.AL(od)+0.01
                Satisficing=[Satisficing Temp_RouteIDs(k)];
            end
        end
        
        Choices=zeros(Assignment.parameters(4),1);
        for i_choices=1:Assignment.parameters(4)
            
            % If there are satisficing routes
            if sum(Satisficing)>0
                
                % Constructing the ranges:
                range=[];
                range=[range, 0];
                for k=1:length(Satisficing)
                    range = [range, k/length(Satisficing)];
                end
                
                % Choosing the alternative
                r = rand();
                m1 = (r >= range);
                m2 = (r <= range);
                m2(1:end-1) = m2(2:end);
                m2(end) = 0;
                
                Choices(i_choices)=Satisficing(find(m1 & m2));
                
                clear range m1 m2
                
                % If there are no satisficing alternatives, the users choose the best alternative available. That is,
                % the minimum utility route.
            else
                [aa,bb] = min(Utility(i_samples,:));
                Choices(i_choices)=Temp_RouteIDs(bb);
            end
            
        end
        
        for i_routes=1:length(Temp_RouteIDs)
            DistributionPerRoute(Temp_RouteIDs(i_routes))=sum(Choices(:)==i_routes)/Assignment.parameters(4);
        end
        
        %                 % RT-DUE
        %             elseif Assignment.Behavior==3
        %
        %                 % Perceived regret.
        %                 Regret=zeros(1,length(Temp_RouteIDs));
        %
        %                 % Calculating the perceived regret:
        %                 Regret=Utility+exp(-Assignment.RTparam.*Utility).*(Utility-min(Utility));
        %
        %                 % Assign users based on the minimization regret behavior:
        %                 for k=1:length(Temp_RouteIDs)
        %                     if (Regret(k)-min(Regret))<0.1 % To avoid numerical fluctuations
        %                         DistributionPerRoute(Temp_RouteIDs(k)) = 1/length(find(Regret - min(Regret) < 0.1));
        %                     else
        %                         DistributionPerRoute(Temp_RouteIDs(k)) = 0;
        %                     end
        %                 end
        
    end
    
    clear Utility k macro_path path_length path_speed i_res Satisficing range r m1 m2 Choices aa bb i_choices Regret
    
    
    
end

ODflowdistrib = DistributionPerRoute;

end


