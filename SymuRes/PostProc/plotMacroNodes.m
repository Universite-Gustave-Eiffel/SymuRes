function plotMacroNodes(Link,Reservoir,ResList,coloringres,MacroNode,NodesList,sizingnodes,Route,RoutesList,coloringroutes,opts)
% plotResMacroNodes(Link,Reservoir,ResList,coloringres,MacroNode,NodesList,sizingnodes,Route,RoutesList,coloringroutes,opts)
% Plot the real network configuration with reservoirs and a set of given
% routes represented by smooth lines. The line thickness represent the
% demand on the route. The routes are shown as sequences of macro nodes.
%
% INPUTS
%---- Link           : Link structure, put [] if no links to plot
%---- Reservoir      : Reservoir structure
%---- ResList        : vector, reservoir IDs
%---- coloringres    : boolean, 1: different colors for the reservoirs
%---- MacroNode      : MacroNode structure
%---- NodesList      : vector, macro nodes IDs
%---- sizingnodes    : boolean, 1: node size depends on the flow transfered
%---- Route          : Route structure
%---- RoutesList     : vector, route IDs
%---- coloringroutes : boolean, 1: different colors for the routes
%---- opts           : options, structure with fields 'fontname', 'fontsize',
%                      'linewidth', 'colormap', 'rescolor', 'textcolor',
%                      'plotlegend', 'plotnumnodes'

NbL = length(Link);
NbR = length(ResList);
NbN = length(NodesList);
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
    FS = 28; % default
end
if isfield(opts,'linewidth')
    LW = opts.linewidth;
else
    LW = 2; % default
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
if isfield(opts,'plotnumnodes')
    plotnumnodes = opts.plotnumnodes;
else
    plotnumnodes = 0; % default
end

% Lines
line0 = {'-', '--', ':', '-.'};
% Line width
minLW = 0.2;
maxLW = 5;
% Marker size
minMS = 5;
maxMS = 30;

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

% Route colors
if coloringroutes == 1
    cmap_routes = arrayextension(cmap0,NbRoutes,'row');
else
    cmap_routes = ones(NbRoutes,1)*rescolor;
end

hold on

% Plot the links
if ~isempty(Link)
    xLinks = zeros(1,2*NbL);
    yLinks = zeros(1,2*NbL);
    for k = 1:NbL
        colorLink = 0.5*[1 1 1];
        plot(Link(k).Points(1,:),Link(k).Points(2,:),'-','Color',colorLink,'LineWidth',LW);
        xLinks(1+2*(k-1)) = Link(k).Points(1,1);
        xLinks(2*k) = Link(k).Points(1,2);
        yLinks(1+2*(k-1)) = Link(k).Points(2,1);
        yLinks(2*k) = Link(k).Points(2,2);
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
mindem = 0.1*maxdem;

arrowL = 0.04*(max(xLinks) - min(xLinks));
hp = zeros(1,NbRoutes);
strleg = cellstr(int2str(zeros(NbRoutes,1)));
i = 1;
for iroute = RoutesList
    if routedem(i) > 0
        sline = '-';
    else
        sline = '--';
    end
    colori = cmap_routes(i,:);
    listx = [];
    listy = [];
    for inode = Route(iroute).NodePath
        xn = MacroNode(inode).Coord(1);
        yn = MacroNode(inode).Coord(2);
        listx = [listx xn];
        listy = [listy yn];
    end
    
    % Smooth the route line
    if length(listx) == 1 % one point: internal trip
        xn = listx(1);
        yn = listy(1);
        xb = mean(Reservoir(r).BorderPoints(1,:));
        yb = mean(Reservoir(r).BorderPoints(2,:));
        d = sqrt((xn - xb)^2 + (yn - yb)^2); % centroid-to-border mean distance
        thmax = 7*pi/4;
        th = 0:0.05:thmax;
        xpath = xn + 0.7*d.*th./thmax.*cos(th);
        ypath = yn + 0.7*d.*th./thmax.*sin(th);
    else
        alpha1 = 0.5; % for way-back turns
        alpha2 = 1.7; % for direct turns
        [xpath, ypath] = smoothroute(listx,listy,50,alpha1,alpha2);
    end
    LWroute = minLW + (routedem(i) - mindem)/(maxdem - mindem)*(maxLW - minLW);
    LWroute = max([LWroute minLW]);
    hp(i) = plot(xpath,ypath,'linestyle',sline,'color',colori,'LineWidth',LWroute);
    strleg{i} = [int2str(iroute) ': [' int2str(Route(iroute).ResPath) ']'];
    singlearrow([xpath(end-1) xpath(end)],[ypath(end-1) ypath(end)],arrowL,'absolute',1,colori,LWroute)
    
    i = i + 1;
end

% Plot the macro nodes
color1 = cmap0(3,:);
color2 = cmap0(1,:);
color3 = cmap0(4,:);
MS1 = 10;
MS2 = 15;
MS3 = 10;
entrynodeslist = [];
exitnodeslist = [];
bordernodeslist = [];
for i = NodesList
    % Modify node coordinates
    %     if strcmp(MacroNode(i).Type,'border')
    %         r1 = MacroNode(i).ResID(1);
    %         r2 = MacroNode(i).ResID(2);
    %         pt11 = Reservoir(r1).Centroid;
    %         pt12 = Reservoir(r2).Centroid;
    %         for ipt = 1:(size(Reservoir(r1).BorderPoints,2)-1)
    %             pt21 = Reservoir(r1).BorderPoints(:,ipt)';
    %             pt22 = Reservoir(r1).BorderPoints(:,ipt+1)';
    %             [isinter, xinter, yinter] = lineintersect(pt11,pt12,pt21,pt22);
    %             if isinter == 1
    %                 MacroNode(i).Coord = [xinter yinter];
    %             end
    %         end
    %     end
    if strcmp(MacroNode(i).Type,'origin') || strcmp(MacroNode(i).Type,'externalentry')
        entrynodeslist = [entrynodeslist i];
    elseif strcmp(MacroNode(i).Type,'destination') || strcmp(MacroNode(i).Type,'externalexit')
        exitnodeslist = [exitnodeslist i];
    else
        bordernodeslist = [bordernodeslist i];
    end
end
if sizingnodes == 0
    for i = exitnodeslist
        plot(MacroNode(i).Coord(1),MacroNode(i).Coord(2),'o','color',color2,'MarkerFaceColor',color2,'markersize',MS2)
    end
    for i = entrynodeslist
        plot(MacroNode(i).Coord(1),MacroNode(i).Coord(2),'o','color',color1,'MarkerFaceColor',color1,'markersize',MS1)
    end
    for i = bordernodeslist
        plot(MacroNode(i).Coord(1),MacroNode(i).Coord(2),'o','color',color3,'MarkerFaceColor',color3,'markersize',MS3)
    end
else
    nodeflow = zeros(1,length(MacroNode));
    for iroute = 1:length(Route)
        if Route(iroute).AssignCoeff > 0
            r = Route(iroute).ResPath(1);
            i_r = Route(iroute).ResRouteIndex(r);
            inode = Route(iroute).NodePath(1); % entry node
            nodeflow(inode) = nodeflow(inode) + mean(Reservoir(r).InflowPerRoute(i_r,:));
            inode = Route(iroute).NodePath(2); % exit node
            nodeflow(inode) = nodeflow(inode) + mean(Reservoir(r).OutflowPerRoute(i_r,:));
            if length(Route(iroute).ResPath) > 1
                k_r = 2;
                for r = Route(iroute).ResPath(2:end)
                    i_r = Route(iroute).ResRouteIndex(r);
                    inode = Route(iroute).NodePath(k_r+1); % exit node
                    nodeflow(inode) = nodeflow(inode) + mean(Reservoir(r).OutflowPerRoute(i_r,:));
                    k_r = k_r + 1;
                end
            end
        end
    end
    maxflow = max(nodeflow(NodesList));
    minflow = 0.1*maxflow;
    for i = exitnodeslist
        MSnode = minMS + (nodeflow(i) - minflow)/(maxflow - minflow)*(maxMS - minMS);
        MSnode = max([MSnode minMS]);
        plot(MacroNode(i).Coord(1),MacroNode(i).Coord(2),'o','color',color2,'MarkerFaceColor',color2,'markersize',MSnode)
    end
    for i = entrynodeslist
        MSnode = minMS + (nodeflow(i) - minflow)/(maxflow - minflow)*(maxMS - minMS);
        MSnode = max([MSnode minMS]);
        plot(MacroNode(i).Coord(1),MacroNode(i).Coord(2),'o','color',color1,'MarkerFaceColor',color1,'markersize',MSnode)
    end
    for i = bordernodeslist
        MSnode = minMS + (nodeflow(i) - minflow)/(maxflow - minflow)*(maxMS - minMS);
        MSnode = max([MSnode minMS]);
        plot(MacroNode(i).Coord(1),MacroNode(i).Coord(2),'o','color',color3,'MarkerFaceColor',color3,'markersize',MSnode)
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
yborder = 0.1; % increasing factor > 0 for the border spacing along y
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

% Plot the macro node numbers
if plotnumnodes == 1
    similnodes = zeros(NbN,NbN);
    for i = NodesList
        for j = NodesList
            dist = sqrt((MacroNode(i).Coord(1) - MacroNode(j).Coord(1))^2 + (MacroNode(i).Coord(2) - MacroNode(j).Coord(2))^2);
            if dist < 0.01*dx %MacroNode(i).ResID(1) == MacroNode(j).ResID(2) && MacroNode(i).ResID(2) == MacroNode(j).ResID(1)
                similnodes(i,j) = 1; % similar nodes if spatially very close
            end
        end
    end
    pairnodeslist = gatherelements(similnodes); % gather nodes that are very close
    for ipair = 1:length(pairnodeslist)
        i = pairnodeslist{ipair}(1);
        text(MacroNode(i).Coord(1),MacroNode(i).Coord(2),['  ' int2str(i)],...
            'color','k','HorizontalAlignment','left','FontName',fontname,'FontSize',FS)
        if length(pairnodeslist{ipair}) > 1
            for i = pairnodeslist{ipair}(2:end) % plot the node ID on the other side for clarity
                text(MacroNode(i).Coord(1),MacroNode(i).Coord(2),[int2str(i) '  '],...
                    'color','k','HorizontalAlignment','right','FontName',fontname,'FontSize',FS)
            end
        end
    end
end

% Plot the legend
if plotlegend == 1
    xleg = xmin + 0.75*(xmax - xmin);
    yleg = ymin + 0.95*(ymax - ymin);
    wleg = 0.05*(xmax - xmin); % legend symbol width
    hleg = 0.06*(ymax - ymin); % height between symbols
    plot([xleg xleg+wleg],yleg*[1 1],'-k','linewidth',minLW)
    text(xleg+wleg,yleg,[' < ' num2str(mindem,2) ' veh/s'],'VerticalAlignment','middle','HorizontalAlignment','left','FontName',fontname,'FontSize',0.8*FS)
    plot([xleg xleg+wleg],(yleg-hleg)*[1 1],'-k','linewidth',maxLW)
    text(xleg+wleg,yleg-hleg,[' > ' num2str(maxdem,2) ' veh/s'],'VerticalAlignment','middle','HorizontalAlignment','left','FontName',fontname,'FontSize',0.8*FS)
    if sizingnodes == 1
        plot(xleg+0.5*wleg,yleg-2*hleg,'o','color','k','MarkerFaceColor','k','markersize',minMS)
        text(xleg+wleg,yleg-2*hleg,[' < ' num2str(minflow,2) ' veh/s'],'VerticalAlignment','middle','HorizontalAlignment','left','FontName',fontname,'FontSize',0.8*FS)
        plot(xleg+0.5*wleg,yleg-3*hleg,'o','color','k','MarkerFaceColor','k','markersize',maxMS)
        text(xleg+wleg,yleg-3*hleg,['  > ' num2str(maxflow,2) ' veh/s'],'VerticalAlignment','middle','HorizontalAlignment','left','FontName',fontname,'FontSize',0.8*FS)
    end
end

axis([xmin xmax ymin ymax])
daspect([1 1 1])
alpha(0.5);

hold off
xlabel('\itx \rm[m]','FontName',fontname,'FontSize',FS)
ylabel('\ity \rm[m]','FontName',fontname,'FontSize',FS)
set(gca,'Position',[0 0 1 1],'FontName',fontname,'FontSize',FS)
set(gca,'visible','off')

end