%% PLOT RESULTS
%--------------------------------------------------------------------------
% DO NOT LAUNCH THIS SCRIPT AT ONCE! Many plots and videos would appear at
% the same time. Go to each section to see what it does.

clear all
clc

addpath('Utilityfunctions/','FDfunctions/')
addpath('MFDsolver/','Assignment/','UserNetworks/','PostProc/')


% GRAPHIC DESIGN
%---------------
% cmap_perso: blue; green; red; yellow; violet; light blue; light violet; light yellow
cmap_perso = [51 51 255; 0 204 51; 204 0 0; 204 153 0; 153 0 102; 51 153 153; 204 102 204; 204 204 102]/255;
cmap_perso2 = [51 51 255; 0 204 51; 204 0 0; 204 153 0; 153 0 102;    51 153 153; 153 204 102; 153 102 102; 204 204 102; 204 102 204]/255;
% cmap_ifsttar: dark blue; blue; blue-green; green
cmap_ifsttar = [0 83 151; 0 118 189; 0 166 164; 87 171 39]/255;
marker_perso = ['o' '+' '*' '^' 'x' 'd' 'v' 's' '<' '>'];
line_perso = {'-', '--', ':', '-.','-', '--', ':', '-.','-', '--', ':', '-.','-', '--', ':', '-.'};
fontname = 'Times New Roman';


%% Load simulation files and launch post-processing
%--------------------------------------------------------------------------

% Result from simulations SC1 and SC2
%--------------------------------------------------------------------------
ScenarioList = {'SC11','SC12','SC21','SC22'};
% SC11: one trip, demand peak, maxdem
% SC12: one trip, demand peak, decrdem
% SC21: one trip, supply reduction, maxdem
% SC22: one trip, supply reduction, decrdem

iresu = 1;
for iscenario = 1:length(ScenarioList)
    % Acc-based
    Results(iresu).Network = 'SingleRes'; % Choice of a network defined by user
    Results(iresu).Solver = 1; % Choice of the solver. 1: accbased / 2: tripbased
    Results(iresu).Name = ScenarioList{iscenario}; % Simulation name
    Results(iresu).Name2 = [ScenarioList{iscenario} ' acc-based']; % Name to print on the graph legends
    
    if Results(iresu).Solver == 1
        load(['UserNetworks/' Results(iresu).Network '/outputs/Outputs_' Results(iresu).Name '_accbased.mat'])
        PostProc_accbased
    elseif Results(iresu).Solver == 2
        load(['UserNetworks/' Results(iresu).Network '/outputs/Outputs_' Results(iresu).Name '_tripbased.mat'])
        PostProc_tripbased
    end
    Results(iresu).Reservoir = Reservoir;
    Results(iresu).Route = Route;
    Results(iresu).Assignment = Assignment;
    Results(iresu).SimulTime = Simulation.Time;
    Results(iresu).ODmacro = ODmacro;
    
    iresu = iresu + 1;
    
    % Trip-based
    Results(iresu).Network = 'SingleRes'; % Choice of a network defined by user
    Results(iresu).Solver = 2; % Choice of the solver. 1: accbased / 2: tripbased
    Results(iresu).Name = ScenarioList{iscenario}; % Simulation name
    Results(iresu).Name2 = [ScenarioList{iscenario} ' trip-based']; % Name to print on the graph legends
    
    if Results(iresu).Solver == 1
        load(['UserNetworks/' Results(iresu).Network '/outputs/Outputs_' Results(iresu).Name '_accbased.mat'])
        PostProc_accbased
    elseif Results(iresu).Solver == 2
        load(['UserNetworks/' Results(iresu).Network '/outputs/Outputs_' Results(iresu).Name '_tripbased.mat'])
        PostProc_tripbased
    end
    Results(iresu).Reservoir = Reservoir;
    Results(iresu).Route = Route;
    Results(iresu).Assignment = Assignment;
    Results(iresu).SimulTime = Simulation.Time;
    Results(iresu).ODmacro = ODmacro;
    
    iresu = iresu + 1;
end



%% Accumulation, inflow and outflow
%--------------------------------------------------------------------------

FS = 16; % font size
FS1 = 16; % title font size
figindex = 'abcdefghijklmnopqrstuvwxyz';

% ResuList = [1 2 3 4]; % list of results: acc-based vs trip-based, SC1
ResuList = [5 6 7 8]; % list of results: acc-based vs trip-based, SC2


% Plot options
ResList = 1; % list of reservoirs
RoutesList = []; % list of routes, put to [] for not plotting route states
TimeRange = [0 Simulation.Duration]; % [s]
AccRange = [0 1000]; % [veh]
FlowRange = [0 2]; % [veh/s]
PlotResTotalVal = 1; % 1: plot reservoir total states / 0: plot route states only
PlotResInternalDyn = 0; % 1: add plots of internal state (debug) / 0: do not add anything (better when plotting several results)
filename = 'test'; % name for printing
Nfig = length(ResList); % number of subfigures
Nplot = length(RoutesList); % number of plots per subfigure
Nresu = length(ResuList); % number of results

% Color and line options
sline = arrayextension(line_perso,Nresu,'column'); % line styles for results
LW = 2*ones(1,Nresu); % line widths for results
% clight = linspace(0,(Nresu-1)/Nresu,Nresu); % color brightness (0 to 1) for results
clight = zeros(1,Nresu); % color brightness (0 to 1) for results
cmap = arrayextension(cmap_perso2,Nplot,'row'); % colors for routes
cmap0 = arrayextension(cmap_perso2,Nresu,'row'); % color for total values (sum over the routes)
% cmap0 = ones(Nresu,1)*[0 0 0]; % color for total values (sum over the routes)

% Figure and subfigure options
Ncol = 1; % number of columns in the figure
Nrow = (floor(Nfig/Ncol) < Nfig/Ncol).*(floor(Nfig/Ncol) + 1) + (floor(Nfig/Ncol) == Nfig/Ncol).*floor(Nfig/Ncol);
[colindex, rowindex] = ind2sub([Ncol Nrow],1:Nfig);
marginleft = 0.1; % figure relative left margin
marginbottom = 0.2; % figure relative bottom margin
intermarginright = 0.15; % subplot relative right margin
intermarginbottom = 0.05; % subplot relative bottom margin
intermargintop = 0.2; % subplot relative top margin
figwidth = 6.5; % whole figure width [inches]
figheight = Nrow*2.5; % whole figure height [inches]

% Accumulation
%--------------------------------------------------------------------------
figure
for ifig = 1:Nfig
    
    subplot(Nrow,Ncol,ifig)
    hp01 = zeros(1,Nresu);
    hp1 = zeros(Nresu,Nplot);
    strleg01 = cellstr(int2str(zeros(Nresu,1)));
    strleg1 = cellstr(int2str(zeros(Nplot,Nresu)));
    hold on
    for iresu = 1:Nresu
        Temp_t = Results(ResuList(iresu)).SimulTime;
        Temp_ntot = Results(ResuList(iresu)).Reservoir(ResList(ifig)).Acc;
        if PlotResTotalVal == 1
            hp01(iresu) = plot(Temp_t,Temp_ntot,...
                'linestyle',sline{iresu},'color',lightencolor(cmap0(iresu,:),clight(iresu)),'linewidth',LW(iresu));
            strleg01{iresu} = [Results(ResuList(iresu)).Name2 ' Total'];
        end
        for iplot = 1:Nplot
            i_r = Results(ResuList(iresu)).Route(RoutesList(iplot)).ResRouteIndex(ResList(ifig));
            if i_r > 0
                Temp_n = Results(ResuList(iresu)).Reservoir(ResList(ifig)).AccPerRoute(i_r,:);
            else
                Temp_n = zeros(1,length(Temp_t));
            end
            hp1(iresu,iplot) = plot(Temp_t,Temp_n,...
                'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
            strleg1{iplot,iresu} = [Results(ResuList(iresu)).Name2 ' Route ' int2str(RoutesList(iplot))];
        end
    end
    hold off
    
    axis([TimeRange AccRange])
    title(['(' figindex(ifig) ') - \itR_{\rm' int2str(ResList(ifig)) '}'],'FontName',fontname,'FontWeight','Bold','FontSize',FS1)
    if ifig + Ncol > Nfig
        xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        if isempty(RoutesList)
            ylabel('accumulation \itn^r(t) \rm[veh]','FontName',fontname,'FontSize',FS)
        else
            ylabel('accumulation \itn_p^r(t) \rm[veh]','FontName',fontname,'FontSize',FS)
        end
    end
    xj = marginleft + (colindex(ifig) - 1)*(1 - marginleft)/Ncol; % left pos
    yi = marginbottom + (Nrow - rowindex(ifig) + intermarginbottom)*(1 - marginbottom)/Nrow; % bottom pos
    w = (1 - intermarginright)*(1 - marginleft)/Ncol; % width
    h = (1 - intermargintop - intermarginbottom)*(1 - marginbottom)/Nrow; % height
    set(gca,'Position',[xj yi w h],'FontName',fontname,'FontSize',FS)
    if ifig == 1
        hp = [];
        strleg = [];
        for iresu = 1:Nresu
            if PlotResTotalVal == 1
                if ~isempty(hp1)
                    hp = [hp hp01(iresu) hp1(iresu,:)];
                    strleg = [strleg; strleg01(iresu); strleg1(:,iresu)];
                else
                    hp = [hp hp01(iresu)];
                    strleg = [strleg; strleg01(iresu)];
                end
            else
                if ~isempty(hp1)
                    hp = [hp hp1(iresu,:)];
                    strleg = [strleg; strleg1(:,iresu)];
                end
            end
        end
        hleg = gridLegend(hp,Nresu,strleg,'Location','North','FontName',fontname,'FontSize',FS);
        legpos = get(hleg,'Position');
        set(hleg,'Position',[xj 0.01 legpos(3) legpos(4)])
    end
end
set(gcf,'Position',[10 10 1000 700])

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
% print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_acc_' filename '.pdf'])

% Inflow
%--------------------------------------------------------------------------
figure
for ifig = 1:Nfig
    
    subplot(Nrow,Ncol,ifig)
    hp01 = zeros(1,Nresu);
    hp1 = zeros(Nresu,Nplot);
    strleg01 = cellstr(int2str(zeros(Nresu,1)));
    strleg1 = cellstr(int2str(zeros(Nplot,Nresu)));
    hold on
    for iresu = 1:Nresu
        Temp_t = Results(ResuList(iresu)).SimulTime;
        Temp_qintot = Results(ResuList(iresu)).Reservoir(ResList(ifig)).Inflow;
        Temp_V = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MeanSpeed;
        Temp_ntot = Results(ResuList(iresu)).Reservoir(ResList(ifig)).Acc;
        Temp_Lavg = Results(ResuList(iresu)).Reservoir(ResList(ifig)).AvgTripLength;
        Temp_Pc = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MaxProd;
        if PlotResTotalVal == 1
            hp01(iresu) = plot(Temp_t,Temp_qintot,...
                'linestyle',sline{iresu},'color',lightencolor(cmap0(iresu,:),clight(iresu)),'linewidth',LW(iresu));
            strleg01{iresu} = [Results(ResuList(iresu)).Name2 ' Total'];
            if PlotResInternalDyn == 1
                plot(Temp_t,Temp_ntot./Temp_Lavg.*Temp_V,...
                    'linestyle','--','color',lightencolor(cmap0(iresu,:),0.4),'linewidth',3);
                plot(Temp_t,Temp_Pc./Temp_Lavg,...
                    'linestyle',':','color',lightencolor(cmap0(iresu,:),0.6),'linewidth',5);
            end
        end
        for iplot = 1:Nplot
            i_r = Results(ResuList(iresu)).Route(RoutesList(iplot)).ResRouteIndex(ResList(ifig));
            if i_r > 0
                Temp_qin = Results(ResuList(iresu)).Reservoir(ResList(ifig)).InflowPerRoute(i_r,:);
                Temp_n = Results(ResuList(iresu)).Reservoir(ResList(ifig)).AccPerRoute(i_r,:);
                Temp_Ltrip = Results(ResuList(iresu)).Reservoir(ResList(ifig)).TripLengthPerRoute(i_r);
            else
                Temp_qin = zeros(1,length(Temp_t));
                Temp_n = zeros(1,length(Temp_t));
                Temp_Ltrip = zeros(1,length(Temp_t));
            end
            hp1(iresu,iplot) = plot(Temp_t,Temp_qin,...
                'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
            strleg1{iplot,iresu} = [Results(ResuList(iresu)).Name2 ' Route ' int2str(RoutesList(iplot))];
            if PlotResInternalDyn == 1
                plot(Temp_t,Temp_n./Temp_Ltrip.*Temp_V,...
                    'linestyle','--','color',lightencolor(cmap(iplot,:),0.4),'linewidth',3);
                plot(Temp_t,Temp_n./Temp_ntot.*Temp_Pc./Temp_Ltrip,...
                    'linestyle',':','color',lightencolor(cmap(iplot,:),0.6),'linewidth',5);
            end
        end
    end
    hold off
    
    axis([TimeRange FlowRange])
    title(['(' figindex(ifig) ') - \itR_{\rm' int2str(ResList(ifig)) '}'],'FontName',fontname,'FontWeight','Bold','FontSize',FS1)
    if ifig + Ncol > Nfig
        xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        if isempty(RoutesList)
            ylabel('inflow \itq_{\rmin}^r(t) \rm[veh/s]','FontName',fontname,'FontSize',FS)
        else
            ylabel('inflow \itq_{\rmin,\itp}^r(t) \rm[veh/s]','FontName',fontname,'FontSize',FS)
        end
    end
    xj = marginleft + (colindex(ifig) - 1)*(1 - marginleft)/Ncol; % left pos
    yi = marginbottom + (Nrow - rowindex(ifig) + intermarginbottom)*(1 - marginbottom)/Nrow; % bottom pos
    w = (1 - intermarginright)*(1 - marginleft)/Ncol; % width
    h = (1 - intermargintop - intermarginbottom)*(1 - marginbottom)/Nrow; % height
    set(gca,'Position',[xj yi w h],'FontName',fontname,'FontSize',FS)
    if ifig == 1
        hp = [];
        strleg = [];
        for iresu = 1:Nresu
            if PlotResTotalVal == 1
                if ~isempty(hp1)
                    hp = [hp hp01(iresu) hp1(iresu,:)];
                    strleg = [strleg; strleg01(iresu); strleg1(:,iresu)];
                else
                    hp = [hp hp01(iresu)];
                    strleg = [strleg; strleg01(iresu)];
                end
            else
                if ~isempty(hp1)
                    hp = [hp hp1(iresu,:)];
                    strleg = [strleg; strleg1(:,iresu)];
                end
            end
        end
        hleg = gridLegend(hp,Nresu,strleg,'Location','North','FontName',fontname,'FontSize',FS);
        legpos = get(hleg,'Position');
        set(hleg,'Position',[xj 0.01 legpos(3) legpos(4)])
    end
end
set(gcf,'Position',[10 10 1000 700])

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
% print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_inflow_' filename '.pdf'])

% Outflow
%--------------------------------------------------------------------------
figure
for ifig = 1:Nfig
    
    subplot(Nrow,Ncol,ifig)
    hp01 = zeros(1,Nresu);
    hp1 = zeros(Nresu,Nplot);
    strleg01 = cellstr(int2str(zeros(Nresu,1)));
    strleg1 = cellstr(int2str(zeros(Nplot,Nresu)));
    hold on
    for iresu = 1:Nresu
        Temp_t = Results(ResuList(iresu)).SimulTime;
        Temp_qouttot = Results(ResuList(iresu)).Reservoir(ResList(ifig)).Outflow;
        Temp_V = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MeanSpeed;
        Temp_ntot = Results(ResuList(iresu)).Reservoir(ResList(ifig)).Acc;
        Temp_Lavg = Results(ResuList(iresu)).Reservoir(ResList(ifig)).AvgTripLength;
        Temp_Pc = Results(ResuList(iresu)).Reservoir(ResList(ifig)).MaxProd;
        if PlotResTotalVal == 1
            hp01(iresu) = plot(Temp_t,Temp_qouttot,...
                'linestyle',sline{iresu},'color',lightencolor(cmap0(iresu,:),clight(iresu)),'linewidth',LW(iresu));
            strleg01{iresu} = [Results(ResuList(iresu)).Name2 ' Total'];
            if PlotResInternalDyn == 1
                plot(Temp_t,Temp_ntot./Temp_Lavg.*Temp_V,...
                    'linestyle','--','color',lightencolor(cmap0(iresu,:),0.4),'linewidth',3);
                plot(Temp_t,Temp_Pc./Temp_Lavg,...
                    'linestyle',':','color',lightencolor(cmap0(iresu,:),0.6),'linewidth',5);
            end
        end
        for iplot = 1:Nplot
            i_r = Results(ResuList(iresu)).Route(RoutesList(iplot)).ResRouteIndex(ResList(ifig));
            if i_r > 0
                Temp_qout = Results(ResuList(iresu)).Reservoir(ResList(ifig)).OutflowPerRoute(i_r,:);
                Temp_n = Results(ResuList(iresu)).Reservoir(ResList(ifig)).AccPerRoute(i_r,:);
                Temp_Ltrip = Results(ResuList(iresu)).Reservoir(ResList(ifig)).TripLengthPerRoute(i_r);
            else
                Temp_qout = zeros(1,length(Temp_t));
                Temp_n = zeros(1,length(Temp_t));
                Temp_Ltrip = zeros(1,length(Temp_t));
            end
            hp1(iresu,iplot) = plot(Temp_t,Temp_qout,...
                'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
            strleg1{iplot,iresu} = [Results(ResuList(iresu)).Name2 ' Route ' int2str(RoutesList(iplot))];
            if PlotResInternalDyn == 1
                plot(Temp_t,Temp_n./Temp_Ltrip.*Temp_V,...
                    'linestyle','--','color',lightencolor(cmap(iplot,:),0.4),'linewidth',3);
                plot(Temp_t,Temp_n./Temp_ntot.*Temp_Pc./Temp_Ltrip,...
                    'linestyle',':','color',lightencolor(cmap(iplot,:),0.6),'linewidth',5);
            end
        end
    end
    hold off
    
    axis([TimeRange FlowRange])
    title(['(' figindex(ifig) ') - \itR_{\rm' int2str(ResList(ifig)) '}'],'FontName',fontname,'FontWeight','Bold','FontSize',FS1)
    if ifig + Ncol > Nfig
        xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        if isempty(RoutesList)
            ylabel('outflow \itq_{\rmout}^r(t) \rm[veh/s]','FontName',fontname,'FontSize',FS)
        else
            ylabel('outflow \itq_{\rmout,\itp}^r(t) \rm[veh/s]','FontName',fontname,'FontSize',FS)
        end
    end
    xj = marginleft + (colindex(ifig) - 1)*(1 - marginleft)/Ncol; % left pos
    yi = marginbottom + (Nrow - rowindex(ifig) + intermarginbottom)*(1 - marginbottom)/Nrow; % bottom pos
    w = (1 - intermarginright)*(1 - marginleft)/Ncol; % width
    h = (1 - intermargintop - intermarginbottom)*(1 - marginbottom)/Nrow; % height
    set(gca,'Position',[xj yi w h],'FontName',fontname,'FontSize',FS)
    if ifig == 1
        hp = [];
        strleg = [];
        for iresu = 1:Nresu
            if PlotResTotalVal == 1
                if ~isempty(hp1)
                    hp = [hp hp01(iresu) hp1(iresu,:)];
                    strleg = [strleg; strleg01(iresu); strleg1(:,iresu)];
                else
                    hp = [hp hp01(iresu)];
                    strleg = [strleg; strleg01(iresu)];
                end
            else
                if ~isempty(hp1)
                    hp = [hp hp1(iresu,:)];
                    strleg = [strleg; strleg1(:,iresu)];
                end
            end
        end
        hleg = gridLegend(hp,Nresu,strleg,'Location','North','FontName',fontname,'FontSize',FS);
        legpos = get(hleg,'Position');
        set(hleg,'Position',[xj 0.01 legpos(3) legpos(4)])
    end
end
set(gcf,'Position',[10 10 1000 700])

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
% print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_outflow_' filename '.pdf'])


%% Total accumulation and inflow/outflow in reservoirs
%--------------------------------------------------------------------------

FS = 16; % font size

% Plot options
ResList = 1; % list of reservoirs
ResuList = [1]; % list of results
TimeRange = [0 Simulation.Duration]; % [s]
AccRange = [0 1000]; % [veh]
FlowRange = [0 2]; % [veh/s]
filename = 'test'; % name for printing
Nplot = length(ResList); % number of plots per figure
Nresu = length(ResuList); % number of results

% Color and line options
sline = arrayextension(line_perso,Nresu,'column'); % line styles for results
LW = 2*ones(1,Nresu); % line widths for results
clight = linspace(0,(Nresu-1)/Nresu,Nresu); % color brightness (0 to 1) for results
cmap = arrayextension(cmap_perso2,Nplot,'row'); % colors for reservoirs

% Figure and subfigure options
figwidth = 6.5; % whole figure width [inches]
figheight = 2.5; % whole figure height [inches]

% Accumulation
%--------------------------------------------------------------------------
figure
strleg = cellstr(int2str(zeros(Nplot*Nresu,1)));
hp1 = zeros(1,Nplot*Nresu);
hold on
i = 1;
for iresu = 1:Nresu
    Temp_t = Results(ResuList(iresu)).SimulTime;
    for iplot = 1:Nplot
        Temp_n = Results(ResuList(iresu)).Reservoir(ResList(iplot)).Acc;
        hp1(i) = plot(Temp_t,Temp_n,...
            'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
        strleg{i} = [Results(ResuList(iresu)).Name2 ' \itR_{\rm' int2str(ResList(iplot)) '}'];
        i = i + 1;
    end
end
hold off
axis([TimeRange AccRange])
xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
ylabel('accumulation \itn^r(t) \rm[veh]','FontName',fontname,'FontSize',FS)
hleg = legend(hp1,strleg);
set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1000 500])

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
% print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_acctot_' filename '.pdf'])

% Inflow
%--------------------------------------------------------------------------
figure
strleg = cellstr(int2str(zeros(Nplot*Nresu,1)));
hp1 = zeros(1,Nplot*Nresu);
hold on
i = 1;
for iresu = 1:Nresu
    Temp_t = Results(ResuList(iresu)).SimulTime;
    for iplot = 1:Nplot
        Temp_qin = Results(ResuList(iresu)).Reservoir(ResList(iplot)).Inflow;
        hp1(i) = plot(Temp_t,Temp_qin,...
            'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
        strleg{i} = [Results(ResuList(iresu)).Name2 ' \itR_{\rm' int2str(ResList(iplot)) '}'];
        i = i + 1;
    end
end
hold off
axis([TimeRange FlowRange])
xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
ylabel('inflow \itq_{\rmin}^r(t) \rm[veh/s]','FontName',fontname,'FontSize',FS)
hleg = legend(hp1,strleg);
set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1000 500])

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
% print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_inflowtot_' filename '.pdf'])

% Outflow
%--------------------------------------------------------------------------
figure
strleg = cellstr(int2str(zeros(Nplot*Nresu,1)));
hp1 = zeros(1,Nplot*Nresu);
hold on
i = 1;
for iresu = 1:Nresu
    Temp_t = Results(ResuList(iresu)).SimulTime;
    for iplot = 1:Nplot
        Temp_qout = Results(ResuList(iresu)).Reservoir(ResList(iplot)).Outflow;
        hp1(i) = plot(Temp_t,Temp_qout,...
            'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
        strleg{i} = [Results(ResuList(iresu)).Name2 ' \itR_{\rm' int2str(ResList(iplot)) '}'];
        i = i + 1;
    end
end
hold off
axis([TimeRange FlowRange])
xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
ylabel('outflow \itq_{\rmout}^r(t) \rm[veh/s]','FontName',fontname,'FontSize',FS)
hleg = legend(hp1,strleg);
set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1000 500])

% set(gcf,'PaperUnits','inches')
% set(gcf,'PaperSize',[figwidth figheight],'PaperPosition',[0 0 figwidth figheight])
% print('-painters','-dpdf',['UserNetworks/' Results(1).Network '/img/' Results(iresu).Network '_outflowtot_' filename '.pdf'])




