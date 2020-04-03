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

% Result from simulation 1
%--------------------------------------------------------------------------
iresu = 1;
Results(iresu).Network = 'SingleRes'; % Choice of a network defined by user
Results(iresu).Solver = 1; % Choice of the solver. 1: accbased / 2: tripbased
Results(iresu).Name = 'SC11'; % Simulation name
Results(iresu).Name2 = 'acc-based'; % Name to print on the graph legends

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


%% Reservoir config and states
%--------------------------------------------------------------------------

% Plot reservoir schematic representation (borders and adjacent connections)
%--------------------------------------------------------------------------
iresu = 1;
figure
opts.fontname = 'Arial';
opts.fontsize = 16;
opts.linewidth = 2;
opts.colormap = cmap_perso2;
plotResBallConfig(Results(iresu).Reservoir,1:NumRes,1,0.5,[1 1],opts)
filename = '';
set(gcf,'PaperUnits','Inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4])
% print('-painters','-dpdf',['UserNetworks/' Results(iresu).Network '/img/' Results(iresu).Network '_resconfig' filename '.pdf'])

% Plot reservoir configuration, real network
%--------------------------------------------------------------------------
iresu = 1;
figure
opts.fontname = 'Arial';
opts.fontsize = 18;
opts.linewidth = 0.8;
opts.colormap = cmap_perso2;
plotResNetConfig([],Results(iresu).Reservoir,1:NumRes,opts)
% filename = '_2';
% set(gcf,'PaperUnits','Inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4])
% print('-painters','-dpdf',['UserNetworks/' Results(iresu).Network '/img/' Results(iresu).Network '_resconfig' filename '.pdf'])

% Plot reservoir state at t
%--------------------------------------------------------------------------
% Plot reservoir state (total accumulation) at t, schematic representation
figure
t0 = 2000;
opts.title = ['\bft = ' num2str(t0) ' s'];
plotResBallAcc(t0,Results(iresu).Reservoir,1:NumRes,Results(iresu).SimulTime,0.5,[1 1],'trafficolor',opts)

% Plot reservoir state (accumulation per route) at t, schematic representation
figure
t0 = 2000;
opts.title = ['\bft = ' num2str(t0) ' s'];
opts.showleg = 1;
plotResBallAccPerRoute(t0,Results(iresu).Reservoir,1:NumRes,Results(iresu).Route,1,Results(iresu).SimulTime,0.5,[1.5 1.5],opts)

% Plot reservoir state (accumulation per route) at t, real network
figure
t0 = 2000;
opts.title = ['\bft = ' num2str(t0) ' s'];
opts.showleg = 1;
plotResNetAccPerRoute(t0,[],Results(iresu).Reservoir,1:NumRes,Results(iresu).Route,1,Results(iresu).SimulTime,opts)

% Plot reservoir state (mean speed) at t, real network
%--------------------------------------------------------------------------
iresu = 1;
figure
t0 = 2000;
SpeedRange = [3 14];
opts.fontname = 'Arial';
opts.fontsize = 16;
opts.linewidth = 1;
opts.title = ['\bfMean speed at t = ' num2str(t0) ' s'];
plotResNetSpeed(t0,[],Results(iresu).Reservoir,1:NumRes,Results(iresu).SimulTime,SpeedRange,'trafficolor',opts)
% filename = '';
% set(gcf,'PaperSize',[11 9],'PaperPosition',[0 0 11 9])
% print('-painters','-dpdf',['UserNetworks/' Results(iresu).Network '/img/' Results(iresu).Network '_meanspeed' filename '.pdf'])

% Plot reservoir total number or demand of routes
%--------------------------------------------------------------------------
figure
opts.legloc = 'East';
plotResRouteDem([],Results(iresu).Reservoir,1:NumRes,Results(iresu).Route,1:NumRoutes,'demand','trafficolor',opts)
% filename = '';
% set(gcf,'PaperSize',[11 9],'PaperPosition',[0 0 11 9])
% print('-painters','-dpdf',['UserNetworks/' Results(iresu).Network '/img/' Results(iresu).Network '_routesdem' filename '.pdf'])

% Plot route paths
%--------------------------------------------------------------------------
iresu = 1;
RoutesList = [];
for iroute = 1:NumRoutes
    if Route(iroute).AssignCoeff > 0 && mean(Results(iresu).Route(iroute).Demand) > 0.0001
        RoutesList = [RoutesList iroute];
    end
end
figure
opts.fontname = 'Arial';
opts.fontsize = 18;
opts.linewidth = 0.8;
opts.plotlegend = 1;
plotRoutes([],Results(iresu).Reservoir,1:NumRes,0,Results(iresu).Route,RoutesList,1,0,opts)
% filename = '';
% set(gcf,'PaperSize',[11 9],'PaperPosition',[0 0 11 9])
% print('-painters','-dpdf',['UserNetworks/' Results(iresu).Network '/img/' Results(iresu).Network '_routes' filename '.pdf'])

% Plot macro nodes
%--------------------------------------------------------------------------
figure
iresu = 1;
opts.fontname = 'Arial';
opts.fontsize = 18;
opts.linewidth = 0.8;
opts.plotlegend = 0;
opts.plotnumnodes = 1;
MacroNodesList = 1:NumMacroNodes;
plotMacroNodes([],Results(iresu).Reservoir,1:NumRes,0,MacroNode,MacroNodesList,0,Results(iresu).Route,[],0,opts)
% filename = '';
% set(gcf,'PaperSize',[11 9],'PaperPosition',[0 0 11 9])
% print('-painters','-dpdf',['UserNetworks/' Results(iresu).Network '/img/' Results(iresu).Network '_macronodes' filename '.pdf'])

% Plot macro nodes with route paths
iresu = 1;
RoutesList = [];
for iroute = 1:NumRoutes
    if Results(iresu).Route(iroute).AssignCoeff > 0 && mean(Results(iresu).Route(iroute).Demand) > 0.0001
        RoutesList = [RoutesList iroute];
    end
end
figure
opts.fontname = 'Arial';
opts.fontsize = 18;
opts.linewidth = 0.8;
opts.plotlegend = 1;
opts.colormap = cmap_perso2;
opts.plotnumnodes = 0;
plotMacroNodes([],Results(iresu).Reservoir,1:NumRes,1,MacroNode,1:NumMacroNodes,1,Results(iresu).Route,RoutesList,0,opts)
% filename = '_routes';
% set(gcf,'PaperUnits','inches','PaperSize',[11 9],'PaperPosition',[0 0 11 9])
% print('-painters','-dpdf',['UserNetworks/' Results(iresu).Network '/img/' Results(iresu).Network '_macronodes' filename '.pdf'])


%% Successive plot of traffic states (video)
%--------------------------------------------------------------------------

% Result to plot
iresu = 1;

% Time frame
tfinal = Simulation.Duration; % final time [s]
Dt_frame = 50; % time step between two frames [s]
Nframes = floor(tfinal/Dt_frame) + 1; % number of frames
t_frame = 0:Dt_frame:tfinal;

% Plot options
opts.fontname = 'Arial';
opts.fontsize = 18;
opts.linewidth = 0.8;
opts.colormap = cmap_perso2;
tref = datenum([2017 1 1 5 00 00]);

figure
for i = 1:Nframes
    clf
    t0 = t_frame(i);
    %opts.title = datestr(tref+t0./(24*3600),'HH:MM');
    opts.title = ['\bft = ' num2str(t0) ' s'];
    opts.showleg = 1;
    %plotResBallAcc(t0,Results(iresu).Reservoir,1:NumRes,Results(iresu).SimulTime,0.5,[1.5 1.5],'trafficolor',opts)
    %plotResBallAccPerRoute(t0,Results(iresu).Reservoir,1:NumRes,Results(iresu).Route,1:NumRoutes,Results(iresu).SimulTime,0.5,[1.5 1.5],opts)
    plotResNetAcc(t0,[],Results(iresu).Reservoir,1:NumRes,Results(iresu).SimulTime,'trafficolor',opts)
    pause(0.001)
    %pause
end



%% Accumulation, inflow and outflow
%--------------------------------------------------------------------------

FS = 16; % font size
FS1 = 16; % title font size
figindex = 'abcdefghijklmnopqrstuvwxyz';

% Plot options
ResList = 1; % list of reservoirs
RoutesList = []; % list of routes, put to [] for not plotting route states
ResuList = [1]; % list of results
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
cmap0 = arrayextension(cmap_perso2,Nplot,'row'); % color for total values (sum over the routes)

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


%% Demand per route
%--------------------------------------------------------------------------

FS = 16; % font size
FS1 = 18; % title font size
figindex = 'abcdefghijklmnopqrstuvwxyz';

% Plot options
RoutesList = 1:NumRoutes; % list of routes
ResuList = [1]; % list of results
TimeRange = [0 Simulation.Duration]; % [s]
FlowRange = [0 2]; % [veh/s]
filename = 'test'; % name for printing
Nplot = length(RoutesList); % number of plots per figure
Nresu = length(ResuList); % number of results

% Color and line options
sline = arrayextension(line_perso,Nresu,'column'); % line styles for results
LW = 2*ones(1,Nresu); % line widths for results
clight = linspace(0,(Nresu-1)/Nresu,Nresu); % color brightness (0 to 1) for results
cmap = arrayextension(cmap_perso2,Nplot,'row'); % colors for routes

% Demand per route
%--------------------------------------------------------------------------
figure
strleg = cellstr(int2str(zeros(Nplot*Nresu,1)));
hp1 = zeros(1,Nplot*Nresu);
hold on
i = 1;
for iresu = 1:Nresu
    for iplot = 1:Nplot
        Temp_tdem = Results(ResuList(iresu)).SimulTime;
        Temp_datadem = Results(ResuList(iresu)).Route(RoutesList(iplot)).Demand;
        Temp_routepath = Results(ResuList(iresu)).Route(RoutesList(iplot)).ResPath;
        hp1(i) = stairs(Temp_tdem,Temp_datadem,...
            'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
        if iplot == 1
            strleg{i} = [Results(ResuList(iresu)).Name2 ' route [' int2str(Temp_routepath) ']'];
        else
            strleg{i} = ['        route [' int2str(Temp_routepath) ']'];
        end
        i = i + 1;
    end
end
hold off
axis([TimeRange FlowRange])
xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
ylabel('demand \rm[veh/s]','FontName',fontname,'FontSize',FS)
hleg = legend(hp1,strleg);
set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1000 500])

% set(gcf,'PaperSize',[30 15],'PaperPosition',[0 0 30 15])
% print('-dpdf',['UserNetworks/' Results(1).Network '/img/demand_comp_' filename '.pdf'])


% Demand per route, subplots
%--------------------------------------------------------------------------
Nfig = length(ResuList); % number of subfigures
Ncol = 2; % number of columns in the figure
Nrow = (floor(Nfig/Ncol) < Nfig/Ncol).*(floor(Nfig/Ncol) + 1) + (floor(Nfig/Ncol) == Nfig/Ncol).*floor(Nfig/Ncol);
Nresu = 1;
figure
for ifig = 1:Nfig
    
    subplot(Nrow,Ncol,ifig)
    strleg = cellstr(int2str(zeros(Nplot*Nresu,1)));
    hp1 = zeros(1,Nplot*Nresu);
    hold on
    i = 1;
    for iresu = 1:Nresu
        for iplot = 1:Nplot
            Temp_tdem = Results(ResuList(iresu)).SimulTime;
            Temp_datadem = Results(ResuList(ifig)).Route(RoutesList(iplot)).Demand;
            Temp_routepath = Results(ResuList(ifig)).Route(RoutesList(iplot)).ResPath;
            hp1(i) = stairs(Temp_tdem,Temp_datadem,...
                'linestyle',sline{iresu},'color',cmap(iplot,:),'linewidth',LW(iresu));
            strleg{i} = ['route [' int2str(Temp_routepath) ']'];
            i = i + 1;
        end
    end
    hold off
    axis([TimeRange FlowRange])
    title(['(' figindex(ifig) ') - ' Results(ifig).Name2],'FontName',fontname,'FontWeight','Bold','FontSize',FS1)
    if ifig + Ncol > Nfig
        xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
    end
    if mod(ifig,Ncol) == 1 || Ncol == 1
        ylabel('demand \rm[veh/s]','FontName',fontname,'FontSize',FS)
    end
    if ifig == 2
        hleg = legend(hp1,strleg);
        set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',0.9*FS)
    end
    set(gca,'FontName',fontname,'FontSize',FS)
end
set(gcf,'Position',[10 10 1000 500])

% set(gcf,'PaperSize',[30 15],'PaperPosition',[0 0 30 15])
% print('-dpdf',['UserNetworks/' Results(1).Network '/img/demand_comp_' filename '_2.pdf'])


%% Demand per OD
%--------------------------------------------------------------------------

FS = 16; % font size
figindex = 'abcdefghijklmnopqrstuvwxyz';

% Plot options
ODList = [1]; % list of OD
ResuList = [1]; % list of results
TimeRange = [0 Simulation.Duration]; % [s]
FlowRange = [0 3]; % [veh/s]
filename = 'test'; % name for printing
Nplot = length(ODList); % number of plots per figure
Nresu = length(ResuList); % number of results

% Color and line options
sline = arrayextension(line_perso,Nresu,'column'); % line styles for results
LW = 2*ones(1,Nresu); % line widths for results
clight = linspace(0,(Nresu-1)/Nresu,Nresu); % color brightness (0 to 1) for results
cmap = arrayextension(cmap_perso2,Nplot,'row'); % colors for routes

figure
strleg = cellstr(int2str(zeros(Nplot*Nresu,1)));
hp1 = zeros(1,Nplot*Nresu);
hold on
i = 1;
for iresu = 1:Nresu
    for iplot = 1:Nplot
        od = ODList(iplot);
        o = ODmacro(od).ResOriginID;
        d = ODmacro(od).ResDestinationID;
        Temp_tdem = [Results(ResuList(iresu)).ODmacro(od).Demand(1).Time  Simulation.Duration];
        Temp_datadem = [Results(ResuList(iresu)).ODmacro(od).Demand(1).Data  Results(ResuList(iresu)).ODmacro(od).Demand(1).Data(end)];
        hp1(i) = stairs(Temp_tdem,Temp_datadem,...
            'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
        strleg{i} = [Results(ResuList(iresu)).Name2 ' {\itR}_{' int2str(o) '} > {\itR}_{' int2str(d) '}'];
        i = i + 1;
    end
end
hold off
axis([TimeRange FlowRange])
xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
ylabel('demand \rm[veh/s]','FontName',fontname,'FontSize',FS)
hleg = legend(hp1,strleg);
set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1000 500])

% set(gcf,'PaperSize',[30 15],'PaperPosition',[0 0 30 15])
% print('-dpdf',['UserNetworks/' Results(1).Network '/img/demandOD_comp_' filename '.pdf'])


%% Travel time per route
%--------------------------------------------------------------------------

FS = 16; % font size

% Plot options
RoutesList = [1 2]; % list of routes
ResuList = [1 2]; % list of results
TimeRange = [0 Simulation.Duration]; % [s]
TTRange = [0 500]; % [s]
PlotFreeflowTT = 0; % 1: add free-flow TT / 0: do not add (better when plotting several results)
filename = 'test'; % name for printing
Nplot = length(RoutesList); % number of plots per figure
Nresu = length(ResuList); % number of results

% Color and line options
sline = arrayextension(line_perso,Nresu,'column'); % line styles for results
LW = 2*ones(1,Nresu); % line widths for results
clight = linspace(0,(Nresu-1)/Nresu,Nresu); % color brightness (0 to 1) for results
cmap = arrayextension(cmap_perso2,Nplot,'row'); % colors for routes

figure
strleg = cellstr(int2str(zeros(Nplot*Nresu,1)));
hp1 = zeros(1,Nplot*Nresu);
hold on
i = 1;
for iresu = 1:Nresu
    for iplot = 1:Nplot
        Temp_t = Results(ResuList(iresu)).SimulTime;
        Temp_TT = Results(ResuList(iresu)).Route(RoutesList(iplot)).TravelTime;
        Temp_TTfree = Results(ResuList(iresu)).Route(RoutesList(iplot)).FreeFlowTravelTime;
        if PlotFreeflowTT == 1
            plot([0 Simulation.Duration],Temp_TTfree*[1 1],...
                'linestyle','--','color',lightencolor(cmap(iplot,:),0.5),'linewidth',4);
        end
        hp1(i) = plot(Temp_t,Temp_TT,...
            'linestyle',sline{iresu},'color',lightencolor(cmap(iplot,:),clight(iresu)),'linewidth',LW(iresu));
        strleg{i} = [Results(ResuList(iresu)).Name2 ' Route ' int2str(RoutesList(iplot)) ': [' int2str(Route(RoutesList(iplot)).ResPath) ']'];
        i = i + 1;
    end
end
hold off
axis([TimeRange TTRange])
xlabel('time \itt \rm[s]','FontName',fontname,'FontSize',FS)
ylabel('travel time \itT(t) \rm[s]','FontName',fontname,'FontSize',FS)
hleg = legend(hp1,strleg);
set(hleg,'Location','best','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1000 500])

% set(gcf,'PaperSize',[30 15],'PaperPosition',[0 0 30 15])
% print('-dpng',['UserNetworks/' Results(1).Network '/img/TT_comp_' filename '.png'])


%% MFD of reservoirs
%--------------------------------------------------------------------------

FS = 16; % font size

% Plot options
ResList = [1 2 3 4]; % list of routes
iresu = 1; % result
AccRange = [0 1000]; % [veh]
ProdRange = [0 3500]; % [veh.m/s]
SpeedRange = [0 15]; % [m/s]
filename = 'test'; % name for printing
Nplot = length(ResList); % number of plots per figure
Nresu = length(ResuList); % number of results

% Color and line options
sline = arrayextension(line_perso,Nplot,'column'); % line styles for reservoirs
LW = linspace(2,5,Nplot); % line widths for reservoirs
cmap = arrayextension(cmap_perso2,Nplot,'row'); % colors for reservoirs


% Production-MFD
%--------------------------------------------------------------------------
figure
strleg = cellstr(int2str(zeros(Nplot,1)));
hp1 = zeros(1,Nplot);
hold on
for iplot = 1:Nplot
    Temp_nj = Results(iresu).Reservoir(ResList(iplot)).MaxAcc;
    Temp_param = Results(iresu).Reservoir(ResList(iplot)).MFDfctParam;
    Temp_n = linspace(0,Temp_nj,100);
    hp1(iplot) = plot(Temp_n,MFDfct(Temp_n,Temp_param),...
        'linestyle',sline{iplot},'color',cmap(iplot,:),'linewidth',LW(iplot));
    strleg{iplot} = ['\itR_{\rm' int2str(ResList(iplot)) '}'];
end
hold off
axis([AccRange ProdRange])
xlabel('accumulation \itn^r \rm[veh]','FontName',fontname,'FontSize',FS)
ylabel('production \itP^r(n^r) \rm[veh.m/s]','FontName',fontname,'FontSize',FS)
hleg = legend(hp1,strleg);
set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1000 500])


% Speed-MFD
%--------------------------------------------------------------------------
figure
strleg = cellstr(int2str(zeros(Nplot,1)));
hp1 = zeros(1,Nplot);
hold on
for iplot = 1:Nplot
    Temp_nj = Results(iresu).Reservoir(ResList(iplot)).MaxAcc;
    Temp_param = Results(iresu).Reservoir(ResList(iplot)).MFDfctParam;
    Temp_n = linspace(0,Temp_nj,100);
    hp1(iplot) = plot(Temp_n,MFDfct(Temp_n,Temp_param)./Temp_n,...
        'linestyle',sline{iplot},'color',cmap(iplot,:),'linewidth',LW(iplot));
    strleg{iplot} = ['\itR_{\rm' int2str(ResList(iplot)) '}'];
end
hold off
axis([AccRange SpeedRange])
xlabel('accumulation \itn^r \rm[veh]','FontName',fontname,'FontSize',FS)
ylabel('mean speed \itV^r(n^r) \rm[m/s]','FontName',fontname,'FontSize',FS)
hleg = legend(hp1,strleg);
set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
set(gca,'FontName',fontname,'FontSize',FS)
set(gcf,'Position',[10 10 1000 500])


%% Video of reservoir states
%--------------------------------------------------------------------------

% Result to plot
iresu = 2;
ResList = 1:NumRes; % list of reservoirs
RoutesList = [1 2]; % list of routes
filename = 'test';

% Time frame
tfinal = Simulation.Duration; % final time [s]
Dt_frame = 200; % time step between two frames [s]
Nframes = floor(tfinal/Dt_frame) + 1; % number of frames
t_frame = 0:Dt_frame:tfinal;

% Plot options
opts.fontname = 'Arial';
opts.fontsize = 18;
opts.linewidth = 0.8;
opts.colormap = cmap_perso2;
tref = datenum([2017 1 1 5 00 00]);

% Make a succession of images
% mkdir(['UserNetworks/' Results(1).Network '/img'],filename)

% Make a video: method 1 (old)
% strfile = ['UserNetworks/' Results(1).Network '/img/' filename '.avi'];
% aviobj = avifile(strfile, 'fps',5);

% Make a video: method 2 (old)
% strfile = ['UserNetworks/' Results(1).Network '/img/' filename '.avi'];
% mov(1:Nframes) = struct('cdata',[], 'colormap',[]);

% Make a video: method 3 (recommended)
% strfile = ['UserNetworks/' Results(1).Network '/img/' filename '.avi'];
% video = VideoWriter(strfile);
% video.FrameRate = 5; % video fps
% open(video)

figure
set(gcf,'Position',[50 50 800 500])
set(gcf,'PaperSize',[30 18],'PaperPosition',[0 0 30 18])
% set(gcf,'PaperSize',[15 11],'PaperPosition',[0 0 15 11])
set(gcf,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');
for i = 1:Nframes
    clf
    t0 = t_frame(i);
    %opts.title = datestr(tref+t0./(24*3600),'HH:MM');
    opts.title = ['\bft = ' num2str(t0) ' s'];
    opts.showleg = 1;
    plotResBallAcc(t0,Results(iresu).Reservoir,ResList,Results(iresu).SimulTime,0.5,[1.5 1.5],'trafficolor',opts)
    %plotResBallAccPerRoute(t0,Results(iresu).Reservoir,ResList,Results(iresu).Route,RoutesList,Results(iresu).SimulTime,0.5,[1.5 1.5],opts)
    %plotResNetAcc(t0,[],Results(iresu).Reservoir,ResList,Results(iresu).SimulTime,'trafficolor',opts)
    %plotResNetAccPerRoute(t0,[],Results(iresu).Reservoir,ResList,Results(iresu).Route,RoutesList,Results(iresu).SimulTime,opts)
    
    % Succession of images
    %print('-dpdf','-painters',['UserNetworks/' Results(1).Network '/img/' filename '/img_' int2str(i) '.pdf'])
    
    % Video method 1 (old)
    %aviobj = addframe(aviobj, getframe(gcf));
    
    % Video method 2 (old)
    %mov(i) = getframe(gcf);
    
    % Video method 3 (recommended)
    %frame = getframe(gcf);
    %writeVideo(video,frame);
    
    pause(0.001)
    %pause
    
end

% Video method 1 (old)
% aviobj = close(aviobj);

% Video method 2 (old)
% movie2avi(mov, strfile, 'compression','Indeo5', 'fps',5);

% Video method 3 (recommended)
% close(video)



%% Video of reservoir states, 2 sub figures
%--------------------------------------------------------------------------

% Result to plot
ResList = 1:NumRes; % list of reservoirs
RoutesList = [1 2]; % list of routes
ResuList = [1 2]; % list of results
filename = 'test';

% Time frame
tfinal = Simulation.Duration; % final time [s]
Dt_frame = 200; % time step between two frames [s]
Nframes = floor(tfinal/Dt_frame) + 1; % number of frames
t_frame = 0:Dt_frame:tfinal;

% Plot options
opts.fontname = 'Arial';
opts.fontsize = 18;
opts.linewidth = 0.8;
opts.colormap = cmap_perso2;
tref = datenum([2017 1 1 5 00 00]);

% Make a succession of images
% mkdir(['UserNetworks/' Results(1).Network '/img'],filename)

% Make a video: method 1 (old)
% strfile = ['UserNetworks/' Results(1).Network '/img/' filename '.avi'];
% aviobj = avifile(strfile, 'fps',5);

% Make a video: method 2 (old)
% strfile = ['UserNetworks/' Results(1).Network '/img/' filename '.avi'];
% mov(1:Nframes) = struct('cdata',[], 'colormap',[]);

% Make a video: method 3 (recommended)
% strfile = ['UserNetworks/' Results(1).Network '/img/' filename '.avi'];
% video = VideoWriter(strfile);
% video.FrameRate = 5; % video fps
% open(video)

Nfig = length(ResuList); % number of subfigures
Ncol = 2; % number of columns in the figure
Nrow = (floor(Nfig/Ncol) < Nfig/Ncol).*(floor(Nfig/Ncol) + 1) + (floor(Nfig/Ncol) == Nfig/Ncol).*floor(Nfig/Ncol);

figure
set(gcf,'Position',[50 50 800 500])
set(gcf,'PaperSize',[30 18],'PaperPosition',[0 0 30 18])
% set(gcf,'PaperSize',[15 11],'PaperPosition',[0 0 15 11])
set(gcf,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');
for i = 1:Nframes
    clf
    t0 = t_frame(i);
    for ifig = 1:Nfig
        subplot(Nrow,Ncol,ifig)
        %opts.title = datestr(tref+t0./(24*3600),'HH:MM');
        opts.title = ['\bft = ' num2str(t0) ' s'];
        opts.showleg = 1;
        %plotResBallAcc(t0,Results(ResuList(ifig)).Reservoir,ResList,Results(ResuList(ifig)).SimulTime,0.5,[1.5 1.5],'trafficolor',opts)
        %plotResBallAccPerRoute(t0,Results(ResuList(ifig)).Reservoir,ResList,Results(ResuList(ifig)).Route,RoutesList,Results(ResuList(ifig)).SimulTime,0.5,[1.5 1.5],opts)
        %plotResNetAcc(t0,[],Results(ResuList(ifig)).Reservoir,ResList,Results(ResuList(ifig)).SimulTime,'trafficolor',opts)
        plotResNetAccPerRoute(t0,[],Results(ResuList(ifig)).Reservoir,ResList,Results(ResuList(ifig)).Route,RoutesList,Results(ResuList(ifig)).SimulTime,opts)
        if ifig == 1
            set(gca,'Position',[0 0.01 0.47 0.95])
        elseif ifig == 2
            set(gca,'Position',[0.55 0.01 0.47 0.95])
        end
    end
    % Succession of images
    %print('-dpdf','-painters',['UserNetworks/' Results(1).Network '/img/' filename '/img_' int2str(i) '.pdf'])
    
    % Video method 1 (old)
    %aviobj = addframe(aviobj, getframe(gcf));
    
    % Video method 2 (old)
    %mov(i) = getframe(gcf);
    
    % Video method 3 (recommended)
    %frame = getframe(gcf);
    %writeVideo(video,frame);
    
    pause(0.001)
    
end

% Video method 1 (old)
% aviobj = close(aviobj);

% Video method 2 (old)
% movie2avi(mov, strfile, 'compression','Indeo5', 'fps',5);

% Video method 3 (recommended)
% close(video)


