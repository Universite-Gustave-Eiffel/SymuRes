function plotLinkReservoirs(Reservoir,ResList,Node,plotnumnodes,Link,plotnumlinks,Entry,Exit,twoways,opts)
% plotLinkReservoirs(Reservoir,ResList,Node,plotnumnodes,Link,plotnumlinks,Entry,Exit,twoways,opts)
% Plot a given list of reservoirs with nodes and links from the real
% network
%
% INPUTS
%---- Reservoir    : Reservoir structure
%---- ResList      : row vector, list of reservoir IDs
%---- Node         : Node structure
%---- plotnumnodes : boolean, plot the node numbers if = 1
%---- Link         : Link structure
%---- plotnumlinks : boolean, plot the link numbers if = 1
%---- Entry        : Entry structure
%---- Exit         : Exit structure
%---- twoways      : boolean, shift the entries and exits numbers if = 1
%                    (for better visibility in a two-way network)
%---- opts         : options, structure with fields 'fontname', 'fontsize',
%                    'linewidth', 'markersize', 'plotaxes', 'colormap'

NbR = length(ResList);
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
    FS = 18; % default
end
if isfield(opts,'linewidth')
    LW = opts.linewidth;
else
    LW = 1; % default
end
if isfield(opts,'markersize')
    MS = opts.markersize;
else
    MS = 5; % default
end
if isfield(opts,'plotaxes')
    plotaxes = opts.plotaxes;
else
    plotaxes = 1; % default
end
if isfield(opts,'colormap')
    cmap = opts.colormap;
else
    cmap = [51 51 255; 0 204 51; 204 0 0; 204 153 0; 153 0 102; 51 153 153; 204 102 204; 204 204 102]/255; % default
end

hold on

ResAdj = cell(1,NbR);
for r = ResList
    ResAdj{r} = Reservoir(r).AdjacentRes;
end
colorRes = vertexcoloring(ResAdj,length(cmap(:,1)));
hp = zeros(1,NbR);
strleg = cellstr(int2str(zeros(NbR,1)));

% Plot the links
if ~isempty(Link)
    
    Npts = 0;
    i = 1;
    for r = ResList
        colorID = colorRes(r);
        colorLink = cmap(colorID,:);
        for k = Reservoir(r).LinksID
            if ~isempty(Link(k).Points)
                hp(i) = plot(Link(k).Points(1,:),Link(k).Points(2,:),'-','Color',colorLink,'LineWidth',LW);
                
                % Plot the link ID
                if plotnumlinks == 1
                    ang = -180 + atan2(((Link(k).Points(2,1))-(Link(k).Points(2,end))),((Link(k).Points(1,1))-(Link(k).Points(1,end))))/pi*180;
                    text(mean(Link(k).Points(1,:)),mean(Link(k).Points(2,:)),int2str(k),'Color',colorLink,...
                        'FontName',fontname,'FontSize',FS,'HorizontalAlignment','center','Rotation',ang,'VerticalAlignment','bottom')
                end
            end
            Npts = Npts + length(Link(k).Points(1,:));
        end
        strleg{i} = ['R ' int2str(r)];
        i = i + 1;
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
    i = 1;
    for r = ResList
        colorID = colorRes(r);
        hp(i) = plot(Reservoir(r).BorderPoints(1,:),Reservoir(r).BorderPoints(2,:),'-','color',colorID,'LineWidth',LW);
        xLinks = [xLinks Reservoir(r).BorderPoints(1,:)];
        yLinks = [yLinks Reservoir(r).BorderPoints(2,:)];
        strleg{i} = ['R ' int2str(r)];
        i = i + 1;
    end
    
end

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
    
    % Plot the nodes
    for r = ResList
        colorID = colorRes(r);
        colorNode = lightencolor(cmap(colorID,:),0.4);
        for i = [Reservoir(r).InternalNodesID Reservoir(r).BoundaryNodesID]
            plot(Node(i).pos0(1),Node(i).pos0(2),'o','color',colorNode,'MarkerFaceColor',colorNode,'markersize',MS)
            
            % Plot the node ID
            if plotnumnodes == 1
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
    end
    
end

% Plot the entries
if ~isempty(Entry)
    
    if twoways == 1
        VertAlign = 'bottom';
    else
        VertAlign = 'middle';
    end
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
if ~isempty(Exit)
    
    if twoways == 1
        VertAlign = 'baseline';
    else
        VertAlign = 'middle';
    end
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

% Plot the reservoir IDs
for r = ResList
    colorID = colorRes(r);
    R0 = 0.03*dx;
    th = 0:0.01:(2*pi);
    fill(Reservoir(r).Centroid(1)+R0*cos(th),Reservoir(r).Centroid(2)+R0*sin(th),'w','EdgeColor',cmap(colorID,:),'FaceAlpha',1,'LineWidth',LW)
    text(Reservoir(r).Centroid(1),Reservoir(r).Centroid(2),int2str(r),'Color',cmap(colorID,:),...
        'HorizontalAlignment','center','FontWeight','bold','FontName',fontname,'FontSize',FS)
end

% Plot size
xborder = 0.05; % increasing factor > 0 for the border spacing along x
yborder = 0.05; % increasing factor > 0 for the border spacing along x
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
% hleg = legend(hp,strleg);
% set(hleg,'Location','EastOutside','FontName',fontname,'FontSize',FS)

hold off

end