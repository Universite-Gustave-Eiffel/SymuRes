function plotResBallAccPerRoute(t,Reservoir,ResList,Route,RoutesList,SimulTime,ResRadius,coordscale,opts)
% plotResBallAccPerRoute(t,Reservoir,ResList,Route,RoutesList,CurrentTime,ResRadius,coordscale,opts)
% Plot the state of reservoirs at time t (accumulation and flow)
% Plot accumulation ratio of each route in the reservoirs
%
% INPUTS
%---- t           : scalar, time [s]
%---- Reservoir   : Reservoir structure
%---- ResList     : vector, reservoir IDs
%---- Route       : Route structure
%---- RoutesList  : vector, route IDs
%---- SimulTime   : vector, simulation time [s]
%---- ResRadius   : scalar, radius of the disk used to symbolize reservoirs
%---- coordscale  : 1-by-2 vector, scale factors along x and y axis [fx fy]
%---- opts        : options, structure with fields 'fontname', 'fontsize',
%                   'colormap', 'textcolor', 'title', 'showleg', 'legloc'

NumRes = length(ResList);
NbRoutes = length(RoutesList);

% timeID = floor(t/TimeStep) + 1; % index of the current time
timeID = findindex(SimulTime,t); % index of the current time

% Options
if isfield(opts,'fontname')
    fontname = opts.fontname;
else
    fontname = 'Arial'; % default
end
if isfield(opts,'fontsize')
    FS = opts.fontsize;
else
    FS = 16; % default
end
if isfield(opts,'colormap')
    cmap0 = opts.colormap;
else
    cmap0 = [51 51 255; 0 204 51; 204 0 0; 204 153 0; 153 0 102; 51 153 153; 204 102 204; 204 204 102]/255; % default
end
if isfield(opts,'textcolor')
    txtcolor = opts.textcolor;
else
    txtcolor = [0.9 0.9 1]; % default
end
if isfield(opts,'title')
    tlabel = opts.title;
else
    tlabel = ['\itt \rm\bf= ' num2str(t) ' s']; % default
end
if isfield(opts,'showleg')
    showleg = opts.showleg;
else
    showleg = 0; % default
end
if isfield(opts,'legloc')
    legloc = opts.legloc;
else
    legloc = 'EastOutside'; % default
end

cmap = cmap0;
while size(cmap,1) < NbRoutes
    cmap = vertcat(cmap,cmap0);
end

flowspac = 0.2; % spacing between flow lines

maxwidth = 30; % flow line max width

xlist = zeros(1,NumRes);
ylist = zeros(1,NumRes);

% Normalization of reservoir coordinates
x0 = Reservoir(ResList(1)).Centroid(1);
y0 = Reservoir(ResList(1)).Centroid(2);
x1 = Reservoir(Reservoir(ResList(1)).AdjacentRes(1)).Centroid(1);
y1 = Reservoir(Reservoir(ResList(1)).AdjacentRes(1)).Centroid(2);
dx0 = max([abs(x1 - x0) abs(y1 - y0)]);

% Define max flow for plotting purpose
listmaxflow = zeros(1,NumRes);
for r = ResList
    if isempty(find(Reservoir(r).AvgTripLength == 0, 1))
        listmaxflow(r) = max(Reservoir(r).MaxProd./Reservoir(r).AvgTripLength);
    end
end
maxflow = max(listmaxflow);

hold on
hf = zeros(1,NbRoutes);
strleg = cellstr(int2str(zeros(NbRoutes,1)));
ileg = 1;
LegList = [];
for r = ResList
    xr = coordscale(1)*(Reservoir(r).Centroid(1) - x0)/dx0;
    yr = coordscale(2)*(Reservoir(r).Centroid(2) - y0)/dx0;
    xlist(r) = xr;
    ylist(r) = yr;
    
    % Plot flow exchanges
    if r < NumRes % to avoid flow line duplication
        for r2 = Reservoir(r).AdjacentRes
            if r2 > r && ismember(r2,ResList) % to avoid flow line duplication
                xr2 = coordscale(1)*(Reservoir(r2).Centroid(1) - x0)/dx0;
                yr2 = coordscale(2)*(Reservoir(r2).Centroid(2) - y0)/dx0;
                ang = atan2(yr2-yr,xr2-xr) + pi/2;
                dx = flowspac*cos(ang);
                dy = flowspac*sin(ang);
                xr_1 = xr + dx;
                yr_1 = yr + dy;
                xr_2 = xr - dx;
                yr_2 = yr - dy;
                xr2_1 = xr2 + dx;
                yr2_1 = yr2 + dy;
                xr2_2 = xr2 - dx;
                yr2_2 = yr2 - dy;
                %maxflow = min([Reservoir(r).MaxProd/Reservoir(r).AvgTripLength(1) Reservoir(r2).MaxProd/Reservoir(r2).AvgTripLength(1)]); % default value
                % Effective flow from Ri to Rj
                flow = sum(Reservoir(r).OutflowPerResPerDest(r2,:,timeID));
                LW = max([flow/maxflow*maxwidth 0.1]);
                %singlearrowtext([xr2_1 xr_1],[yr2_1 yr_1],'',0.8,0.1,'k',LW,fontname,FS,'center','middle')
                plot([xr_1 xr2_1],[yr_1 yr2_1],'-','color','k','linewidth',LW)
                text(1/3*xr_1+2/3*xr2_1,1/3*yr_1+2/3*yr2_1,num2str(flow,2),'Rotation',ang*180/pi,...
                    'HorizontalAlignment','center','color','k','BackgroundColor','w','FontName',fontname,'FontSize',0.5*FS)
                % Effective flow from Rj to Ri
                flow = sum(Reservoir(r2).OutflowPerResPerDest(r,:,timeID));
                LW = max([flow/maxflow*maxwidth 0.1]);
                plot([xr_2 xr2_2],[yr_2 yr2_2],'-','color','k','linewidth',LW)
                text(2/3*xr_2+1/3*xr2_2,2/3*yr_2+1/3*yr2_2,num2str(flow,2),'Rotation',ang*180/pi,...
                    'HorizontalAlignment','center','color','k','BackgroundColor','w','FontName',fontname,'FontSize',0.5*FS)
            end
        end
    end
    
    % Plot reservoir disk
    th = 0:0.01:(2*pi);
    x = xr + ResRadius*cos(th);
    y = yr + ResRadius*sin(th);
    fill(x,y,0.6*[1 1 1],'EdgeColor','none')
    
    % Plot accumulation evolution
    angstart = 0;
    k_r = 1;
    for iroute = RoutesList
        i_r = find(Reservoir(r).RoutesID == iroute);
        if ~isempty(i_r)
            accratio = Reservoir(r).AccPerRoute(i_r,timeID)/Reservoir(r).MaxAcc;
            angend = angstart + accratio*2*pi;
            th = angstart:0.01:angend;
            x = xr + ResRadius*cos(th);
            y = yr + ResRadius*sin(th);
            angstart = angend;
            if ~ismember(iroute,LegList) % add to the legend
                hf(ileg) = fill([xr x],[yr y],cmap(k_r,:),'EdgeColor','none');
                strleg{ileg} = ['[' int2str(Route(iroute).ResPath) ']'];
                LegList = [LegList iroute];
                ileg = ileg + 1;
            else
                fill([xr x],[yr y],cmap(k_r,:),'EdgeColor','none');
            end
        end
        k_r = k_r + 1;
    end
    
    text(xr,yr,{['{\itR}_{' int2str(r) '}']; int2str(round(Reservoir(r).Acc(timeID)))},...
        'HorizontalAlignment','center','color',txtcolor,'FontName',fontname,'FontWeight','Bold','FontSize',FS)
end

% Plot size
xborder = 0.05; % increasing factor > 0 for the border spacing along x
yborder = 0.1; % increasing factor > 0 for the border spacing along x
if max(xlist) == min(xlist)
    dx = max(ylist) - min(ylist);
else
    dx = max(xlist) - min(xlist);
end
if max(ylist) == min(ylist)
    dy = max(xlist) - min(xlist);
else
    dy = max(ylist) - min(ylist);
end
xmin = min(xlist) - ResRadius - xborder*dx;
xmax = max(xlist) + ResRadius + xborder*dx;
ymin = min(ylist) - ResRadius - yborder*dy;
ymax = max(ylist) + ResRadius + yborder*dy;

% text((min(xlist)+max(xlist))/2,max(ylist)+2*ResRadius,{['\itt \rm\bf= ' num2str(t) ' s']; '\rmAccumulation [veh]'},...
%     'HorizontalAlignment','center','FontName',fontname,'FontWeight','Bold','FontSize',FS)
text((xmin+xmax)/2,ymax,tlabel,...
    'HorizontalAlignment','center','VerticalAlignment','top','FontName',fontname,'FontSize',FS)
hold off
axis([xmin xmax ymin ymax])
daspect([1 1 1])
set(gca,'Position',[0 0 1 1],'FontName',fontname,'FontSize',FS)
set(gca,'visible','off')
if showleg == 1
    hleg = legend(hf,strleg);
    set(hleg,'Location',legloc,'FontName',fontname,'FontSize',FS)
end

end

