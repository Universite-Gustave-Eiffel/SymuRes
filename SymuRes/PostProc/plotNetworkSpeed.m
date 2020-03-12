function plotNetworkSpeed(Node,Link,Entry,Exit,opts)
% plotNetwork(Node,Link,Entry,Exit,opts)
% Plot the real network configuration with link speeds
%
% INPUTS
%---- Node  : Node structure
%---- Link  : Link structure
%---- Entry : Entry structure
%---- Exit  : Exit structure
%---- opts  : options, structure with fields 'fontname', 'fontsize',
%             'linewidth', 'markersize', 'plotnumnodes', 'plotnumlinks',
%             'twoways', 'plotaxes', 'colormap', 'speedtime', 'speedrange'

NbN = length(Node);
NbL = length(Link);
NbE = length(Entry);
NbS = length(Exit);

% Options
if isfield(opts,'fontname')
    fontname = opts.fontname;
else
    fontname = 'Arial'; % default
end
if isfield(opts,'fontsize')
    FS = opts.fontsize;
else
    FS = 12; % default
end
if isfield(opts,'linewidth')
    LW = opts.linewidth;
else
    LW = 2; % default
end
if isfield(opts,'markersize')
    MS = opts.markersize;
else
    MS = 5; % default
end
if isfield(opts,'plotnumnodes')
    plotnumnodes = opts.plotnumnodes;
else
    plotnumnodes = 0; % default
end
if isfield(opts,'plotnumlinks')
    plotnumlinks = opts.plotnumlinks;
else
    plotnumlinks = 0; % default
end
if isfield(opts,'twoways')
    twoways = opts.twoways;
else
    twoways = 0; % default
end
if isfield(opts,'plotaxes')
    plotaxes = opts.plotaxes;
else
    plotaxes = 1; % default
end
if isfield(opts,'colormap')
    cmap = opts.colormap;
else
    cmap = [0 0 0; 0.8 0 0; 0 0.25 1; 1 0.25 0]; % default
end
if isfield(opts,'speedtime')
    speedtime = opts.speedtime;
else
    speedtime = 0; % default, do not show link speed
end
if isfield(opts,'speedrange')
    speedrange = opts.speedrange;
else
    speedrange = [0 50/3.6]; % default [m/s]
end

hold on

% Plot the links
if ~isempty(Link)
    Npts = 0;
    for k = 1:NbL
        if speedtime > 0
            speedratio = (Link(k).MeanSpeed(speedtime) - speedrange(1))/(speedrange(2) - speedrange(1));
            colormap(cmap)
            nbColor = size(cmap,1);
            indcolor = min([max([floor(speedratio*nbColor) 1]) nbColor]);
            colorLink = cmap(indcolor,:);
        else
            colorLink = cmap(1,:);
        end
        plot(Link(k).Points(1,:),Link(k).Points(2,:),'-','Color',colorLink,'LineWidth',LW);
        ang = -180 + atan2(((Link(k).Points(2,1))-(Link(k).Points(2,end))),((Link(k).Points(1,1))-(Link(k).Points(1,end))))/pi*180;
        if plotnumlinks == 1
            text(mean(Link(k).Points(1,:)),mean(Link(k).Points(2,:)),int2str(k),'Color',colorLink,...
                'FontName',fontname,'FontSize',FS,'HorizontalAlignment','center','Rotation',ang,'VerticalAlignment','bottom')
        end
        Npts = Npts + length(Link(k).Points(1,:));
    end
    xLinks = zeros(1,Npts);
    yLinks = zeros(1,Npts);
    i = 1;
    for k = 1:NbL
        Nk = length(Link(k).Points(1,:));
        xLinks(i:(i+Nk-1)) = Link(k).Points(1,:);
        yLinks(i:(i+Nk-1)) = Link(k).Points(2,:);
        i = i + Nk;
    end
end

% Add a colorbar for link speed values
if speedtime > 0
    hcb = colorbar('East');
    hcb.Label.String = 'Mean speed [km/h]';
    hcb.Label.FontSize = FS;
    hcb.Label.FontName = fontname;
    caxis(3.6*speedrange);
end

% Plot the nodes
if ~isempty(Node)
    % Append coordinates and position to Node structure
    for i = 1:NbN
        Node(i).Coord = [];
    end
    for k=1:NbL
        for j=1:length(Link(k).NodeUpID)
            Node(Link(k).NodeUpID(j)).Coord = [Node(Link(k).NodeUpID(j)).Coord Link(k).Points(:,1)];
        end
        for j=1:length(Link(k).NodeDownID)
            Node(Link(k).NodeDownID(j)).Coord = [Node(Link(k).NodeDownID(j)).Coord Link(k).Points(:,end)];
        end
    end
    for i = 1:NbN
        if(size(Node(i).Coord,2)>1)
            Node(i).pos0 = mean(Node(i).Coord,2);
        else
            Node(i).pos0 = Node(i).Coord;
        end
    end
    
    xLinks = zeros(1,NbN);
    yLinks = zeros(1,NbN);
    for i = 1:NbN
        colorNode = cmap(2,:);
        plot(Node(i).pos0(1),Node(i).pos0(2),'o','color',colorNode,'MarkerFaceColor',colorNode,'markersize',MS)
        xLinks(i) = Node(i).pos0(1);
        yLinks(i) = Node(i).pos0(2);
    end
end
if plotnumnodes == 1
    for i = 1:NbN
        
        if Node(i).Type == 1 || Node(i).Type == 2
            if Node(i).Type == 1 % entry node
                idL = Node(i).OutgoingLinksID;
                ang = 180 + atan2(((Link(idL).Points(2,1))-(Link(idL).Points(2,end))),((Link(idL).Points(1,1))-(Link(idL).Points(1,end))))/pi*180;
            elseif Node(i).Type == 2 % exit node
                idL = Node(i).IncomingLinksID;
                ang = 180 + atan2(((Link(idL).Points(2,1))-(Link(idL).Points(2,end))),((Link(idL).Points(1,1))-(Link(idL).Points(1,end))))/pi*180;
            end
            if (0 <= ang && ang < 45) || (315 <= ang && ang <= 360)
                VertAlign = 'bottom';
                HoriAlign = 'center';
            elseif 45 <= ang && ang < 135
                VertAlign = 'middle';
                HoriAlign = 'right';
            elseif 135 <= ang && ang < 225
                VertAlign = 'top';
                HoriAlign = 'center';
            elseif 225 <= ang && ang < 315
                VertAlign = 'middle';
                HoriAlign = 'left';
            end
            
        else % internal node
            VertAlign = 'middle';
            HoriAlign = 'center';
        end
        
        if twoways == 0
            VertAlign = 'middle';
            HoriAlign = 'center';
        end
        text(Node(i).pos0(1),Node(i).pos0(2),int2str(i),'Color',[0 0.5 0],'FontSize',FS,...
            'FontName',fontname,'HorizontalAlignment',HoriAlign,'VerticalAlignment',VertAlign,'FontWeight','bold','BackgroundColor',[1 1 1])
    end
end

% Plot the entries
if twoways == 1
    VertAlign = 'bottom';
else
    VertAlign = 'middle';
end
if ~isempty(Entry)
    for k = 1:NbE
        colorEntry = cmap(3,:);
        xt = Node(Entry(k).NodeID).pos0(1);
        yt = Node(Entry(k).NodeID).pos0(2);
        idL = Node(Entry(k).NodeID).OutgoingLinksID;
        ang = 180 + atan2(((Link(idL).Points(2,1))-(Link(idL).Points(2,end))),((Link(idL).Points(1,1))-(Link(idL).Points(1,end))))/pi*180;
        text(xt,yt,['O' int2str(k) '  \rightarrow  '],'FontWeight','bold','Color',colorEntry,'FontSize',FS,...
            'FontName',fontname,'Rotation',ang,'HorizontalAlignment','right','VerticalAlignment',VertAlign);
    end
end

% Plot the exits
if twoways == 1
    VertAlign = 'baseline';
else
    VertAlign = 'middle';
end
if ~isempty(Exit)
    for k=1:NbS
        colorExit = cmap(4,:);
        xt = Node(Exit(k).NodeID).pos0(1);
        yt = Node(Exit(k).NodeID).pos0(2);
        idL = Node(Exit(k).NodeID).IncomingLinksID;
        ang = 180 + atan2(((Link(idL).Points(2,1))-(Link(idL).Points(2,end))),((Link(idL).Points(1,1))-(Link(idL).Points(1,end))))/pi*180;
        text(xt,yt,['  \rightarrow  D' int2str(k)],'FontWeight','bold','Color',colorExit,'FontSize',FS,...
            'FontName',fontname,'Rotation',ang,'HorizontalAlignment','left','VerticalAlignment',VertAlign);
    end
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

if plotaxes == 1
    xlabel('\itx \rm[m]','FontName',fontname,'FontSize',FS)
    ylabel('\ity \rm[m]','FontName',fontname,'FontSize',FS)
    set(gca,'FontName',fontname,'FontSize',FS)
else
    set(gca,'visible','off','Position',[0 0 1 1])
end
% hold off

end