function plotResNetConfig(Link,Reservoir,ResList,opts)
% plotResNetConfig(Link,Reservoir,ResList,opts)
% Plot the reservoir configuration with links and/or shape borders
%
% INPUTS
%---- Link      : Link structure, put [] if no links to plot
%---- Reservoir : Reservoir structure
%---- ResList   : vector, list of the reservoir IDs
%---- opts      : options, structure with fields 'fontname', 'fontsize',
%                 'linewidth', 'colormap', 'textcolor'

NbL = length(Link);
NbR = length(ResList);

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
    cmap = opts.colormap;
else
    cmap = [51 51 255; 0 204 51; 204 0 0; 204 153 0; 153 0 102; 51 153 153; 204 102 204; 204 204 102]/255; % default
end
if isfield(opts,'textcolor')
    txtcolor = opts.textcolor;
else
    txtcolor = [0.9 0.9 1]; % default
end

% Reservoir colors
ResAdj = cell(1,NbR);
for r = ResList
    ResAdj{r} = intersect(Reservoir(r).AdjacentRes,ResList);
end
colorRes = vertexcoloring(ResAdj,length(cmap(:,1)));

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
    for k = 1:NbL
        Nk = length(Link(k).Points(1,:));
        xLinks(i:(i+Nk-1)) = Link(k).Points(1,:);
        yLinks(i:(i+Nk-1)) = Link(k).Points(2,:);
        i = i + Nk;
    end
else
    % Case when the reservoir borders are defined but not the link network
    xLinks = [];
    yLinks = [];
    for r = ResList
        xLinks = [xLinks Reservoir(r).BorderPoints(1,:)];
        yLinks = [yLinks Reservoir(r).BorderPoints(2,:)];
    end
end

% Plot the reservoirs
i = 1;
for r = ResList
    colori = cmap(colorRes(i),:);
    fill(Reservoir(r).BorderPoints(1,:),Reservoir(r).BorderPoints(2,:),colori,'EdgeColor','none')
    plot(Reservoir(r).BorderPoints(1,:),Reservoir(r).BorderPoints(2,:),'-','color',colori,'LineWidth',LW);
    i = i + 1;
end
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

xlabel('\itx \rm[m]','FontName',fontname,'FontSize',FS)
ylabel('\ity \rm[m]','FontName',fontname,'FontSize',FS)
set(gca,'Position',[0 0 1 1],'FontName',fontname,'FontSize',FS)
set(gca,'visible','off')
hold off

end