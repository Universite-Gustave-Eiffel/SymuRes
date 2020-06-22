function plotRoutes(Link,Reservoir,ResList,coloringres,MacroNode,NodesList,coloringnodes,Route,RoutesList,coloringroutes,sizingroutes,opts)
% plotRoutes(Link,Reservoir,ResList,coloringres,MacroNode,NodesList,coloringnodes,Route,RoutesList,coloringroutes,sizingroutes,opts)
% Plot the real network configuration with reservoirs and a set of given
% routes represented by smooth lines.
%
% INPUTS
%---- Link           : Link structure, put [] if no links to plot
%---- Reservoir      : Reservoir structure
%---- ResList        : vector, reservoir IDs
%---- coloringres    : boolean, 1: different colors for the reservoirs
%---- MacroNode      : MacroNode structure
%---- NodesList      : vector, macro nodes IDs
%---- coloringnodes  : boolean, 1: different colors by node types
%---- Route          : Route structure
%---- RoutesList     : vector, route IDs
%---- coloringroutes : boolean, 1: different colors for the routes
%---- sizingroutes   : boolean, 1: route thickness depends on route mean demand
%---- opts           : options, structure with fields 'fontname', 'fontsize',
%                      'linewidth', 'markersize', 'colormap', 'rescolor', 'textcolor',
%                      'plotlegend', 'exactsmooth', 'nodepath'

NbR = length(ResList);
NbRoutes = length(RoutesList);

% Options
if isfield(opts,'fontname')
    fontname = opts.fontname;
else
    fontname = 'Arial'; % default
end
if isfield(opts,'fontsize')
    FS = opts.fontsize;
else
    FS = 24; % default
end
if isfield(opts,'linewidth')
    LW = opts.linewidth;
else
    LW = 2; % default
end
if isfield(opts,'markersize')
    MS = opts.markersize;
else
    MS = 10; % default
end
if isfield(opts,'colormap')
    cmap0 = opts.colormap;
else
    cmap0 = [51 51 255; 0 204 51; 204 0 0; 204 153 0; 153 0 102; 51 153 153; 204 102 204; 204 204 102]/255; % default
end
if isfield(opts,'rescolor')
    rescolor = opts.rescolor;
else
    rescolor = [0.1 0.1 0]; % default
end
if isfield(opts,'textcolor')
    txtcolor = opts.textcolor;
else
    txtcolor = [0.9 0.9 1]; % default
end
if isfield(opts,'plotlegend')
    plotlegend = opts.plotlegend;
else
    plotlegend = 0; % default
end
if isfield(opts,'exactsmooth')
    exactsmooth = opts.exactsmooth;
else
    exactsmooth = 1; % default
end
if isfield(opts,'nodepath')
    plotnodepath = opts.nodepath;
else
    plotnodepath = 1; % default
end

% Lines
line0 = {'-', '--', ':', '-.'};
% Line width
minLW = 0.2;
maxLW = 5;

% Reservoir colors
if coloringres == 1
    ResAdj = cell(1,NbR);
    for r = ResList
        ResAdj{r} = intersect(Reservoir(r).AdjacentRes,ResList);
    end
    colorIDlist = vertexcoloring(ResAdj,length(cmap0(:,1)));
    cmap_res = cmap0(colorIDlist,:);
else
    cmap_res = ones(NbR,1)*rescolor;
end

% Macro node colors
if coloringnodes == 1
    color1 = cmap0(3,:);
    color2 = cmap0(1,:);
    color3 = cmap0(4,:);
else
    color1 = [0 0 0];
    color2 = [0 0 0];
    color3 = [0 0 0];
end

% Route colors and lines
if coloringroutes == 1
    cmap_routes = arrayextension(cmap0,NbRoutes,'row');
    line_routes = arrayextension(line0,NbRoutes,'column');
else
    cmap_routes = ones(NbRoutes,1)*rescolor;
    line_routes = cell(1,NbRoutes);
    for i = 1:NbRoutes
        line_routes{i} = '-';
    end
end


hold on

% Plot the links
if ~isempty(Link)
    Npts = 0;
    for r = ResList
        for k = Reservoir(r).LinksID
            colorLink = 0.5*[1 1 1];
            plot(Link(k).Points(1,:),Link(k).Points(2,:),'-','Color',colorLink,'LineWidth',LW);
            Npts = Npts + length(Link(k).Points(1,:));
        end
    end
    xLinks = zeros(1,Npts);
    yLinks = zeros(1,Npts);
    i = 1;
    for r = ResList
        for k = Reservoir(r).LinksID
            Nk = length(Link(k).Points(1,:));
            xLinks(i:(i+Nk-1)) = Link(k).Points(1,:);
            yLinks(i:(i+Nk-1)) = Link(k).Points(2,:);
            i = i + Nk;
        end
    end
else
    % Case when the reservoir borders are defined but not the link network
    xLinks = [];
    yLinks = [];
    for r = ResList
        if ~isempty(Reservoir(r).BorderPoints)
            xLinks = [xLinks Reservoir(r).BorderPoints(1,:)];
            yLinks = [yLinks Reservoir(r).BorderPoints(2,:)];
        end
    end
end

% Plot the reservoirs
i = 1;
for r = ResList
    if ~isempty(Reservoir(r).BorderPoints)
        colori = cmap_res(i,:);
        fill(Reservoir(r).BorderPoints(1,:),Reservoir(r).BorderPoints(2,:),colori,'EdgeColor','none')
        plot(Reservoir(r).BorderPoints(1,:),Reservoir(r).BorderPoints(2,:),'-','color',colori,'LineWidth',LW);
    end
    i = i + 1;
end

% Plot the routes
i = 1;
routedem = zeros(1,NbRoutes);
for iroute = RoutesList
    routedem(i) = mean(Route(iroute).Demand);
    i = i + 1;
end
maxdem = max(routedem);

arrowL = 0.04*(max(xLinks) - min(xLinks));
hp = zeros(1,NbRoutes);
strleg = cellstr(int2str(zeros(NbRoutes,1)));
i = 1;
for iroute = RoutesList
    if routedem(i) > 0 || sizingroutes == 0
        colori = cmap_routes(i,:);
        listx = [];
        listy = [];
        if plotnodepath == 1
            for m = Route(iroute).NodePath
                xm = MacroNode(m).Coord(1);
                ym = MacroNode(m).Coord(2);
                listx = [listx xm];
                listy = [listy ym];
            end
        else
            for r = Route(iroute).ResPath
                xr = Reservoir(r).Centroid(1);
                yr = Reservoir(r).Centroid(2);
                listx = [listx xr];
                listy = [listy yr];
            end
        end
        
        % Smooth the route line
        if length(listx) == 1 % one point: internal trip
            xr = listx(1);
            yr = listy(1);
            xb = mean(Reservoir(r).BorderPoints(1,:));
            yb = mean(Reservoir(r).BorderPoints(2,:));
            d = sqrt((xr - xb)^2 + (yr - yb)^2); % centroid-to-border mean distance
            thmax = 7*pi/4;
            th = 0:0.05:thmax;
            xpath = xr + 0.7*d.*th./thmax.*cos(th);
            ypath = yr + 0.7*d.*th./thmax.*sin(th);
        else
            if exactsmooth == 0
                % The smoothed route does not necessarily connect all the points
                alpha1 = 0.5; % for way-back turns
                alpha2 = 1.7; % for direct turns
                [xpath, ypath] = smoothroute(listx,listy,50,alpha1,alpha2);
            else
                % The smoothed route connects all the points (exact interpolation)
                tension = 0.5; % smooth coeff
                [xpath, ypath] = smoothroute2(listx,listy,50,tension);
            end
        end
        if sizingroutes == 1
            LWroute = minLW + routedem(i)/maxdem*(maxLW - minLW);
        else
            LWroute = 2*LW;
        end
        hp(i) = plot(xpath,ypath,'linestyle',line_routes{i},'color',colori,'LineWidth',LWroute);
        strleg{i} = [int2str(iroute) ': [' int2str(Route(iroute).ResPath) ']'];
        singlearrow([xpath(end-1) xpath(end)],[ypath(end-1) ypath(end)],arrowL,'absolute',1,colori,LWroute)
    end
    
    i = i + 1;
end

% Plot the macro nodes
MS1 = MS;
MS2 = 1.5*MS;
MS3 = MS;
for i = NodesList
    if strcmp(MacroNode(i).Type,'origin') || strcmp(MacroNode(i).Type,'externalentry')
        plot(MacroNode(i).Coord(1),MacroNode(i).Coord(2),'o','color',color2,'MarkerFaceColor',color2,'markersize',MS2)
    elseif strcmp(MacroNode(i).Type,'destination') || strcmp(MacroNode(i).Type,'externalexit')
        plot(MacroNode(i).Coord(1),MacroNode(i).Coord(2),'o','color',color1,'MarkerFaceColor',color1,'markersize',MS1)
    else
        plot(MacroNode(i).Coord(1),MacroNode(i).Coord(2),'o','color',color3,'MarkerFaceColor',color3,'markersize',MS3)
    end
end

% Plot the reservoir numbers
for r = ResList
    xr = Reservoir(r).Centroid(1);
    yr = Reservoir(r).Centroid(2);
    text(xr,yr,['{\itR}_{' int2str(r) '}'],...
        'HorizontalAlignment','center','color',txtcolor,'FontName',fontname,'FontWeight','Bold','FontSize',FS)
end

% Plot size
xborder = 0.1; % increasing factor > 0 for the border spacing along x
yborder = 0.1; % increasing factor > 0 for the border spacing along x
if max(xLinks) == min(xLinks)
    dx = max(yLinks) - min(yLinks);
else
    dx = max(xLinks) - min(xLinks);
end
if max(yLinks) == min(yLinks)
    dy = max(xLinks) - min(xLinks);
else
    dy = max(yLinks) - min(yLinks);
end
xmin = min(xLinks) - xborder*dx;
xmax = max(xLinks) + xborder*dx;
ymin = min(yLinks) - yborder*dy;
ymax = max(yLinks) + yborder*dy;

axis([xmin xmax ymin ymax])
daspect([1 1 1])
alpha(0.5);

% text((xmin+xmax)/2,ymax-yborder*dy/2,'Route paths',...
%     'HorizontalAlignment','center','FontName',fontname,'FontWeight','Bold','FontSize',FS)
hold off
xlabel('\itx \rm[m]','FontName',fontname,'FontSize',FS)
ylabel('\ity \rm[m]','FontName',fontname,'FontSize',FS)
set(gca,'Position',[0 0 1 1],'FontName',fontname,'FontSize',FS)
set(gca,'visible','off')
if plotlegend == 1
    hleg = legend(hp,strleg);
    set(hleg,'Location','NorthEast','FontName',fontname,'FontSize',FS)
end

end