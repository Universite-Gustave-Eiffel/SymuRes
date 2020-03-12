function plotResBallConfig(Reservoir,ResList,plotborders,ResRadius,coordscale,opts)
% plotResBallConfig(Reservoir,ResList,plotborders,ResRadius,coordscale,opts)
% Plot the reservoir configuration (ball shapes and connections)
%
% INPUTS
%---- Reservoir   : Reservoir structure
%---- ResList     : vector, reservoir IDs
%---- plotborders : boolean, plot the reservoir borders if = 1
%---- ResRadius   : scalar, radius of the disk used to symbolize reservoirs
%---- coordscale  : 1-by-2 vector, scale factors along x and y axis [fx fy]
%---- opts        : options, structure with fields 'fontname', 'fontsize',
%                   'linewidth', 'colormap', 'textcolor'

NumRes = length(ResList);

% Options
if isfield(opts,'fontname')
    fontname = opts.fontname;
else
    fontname = 'Arial'; % default
end
if isfield(opts,'fontsize')
    FS = opts.fontsize;
else
    FS = 20; % default
end
if isfield(opts,'linewidth')
    LW = opts.linewidth;
else
    LW = 1.5; % default
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
ResAdj = cell(1,NumRes);
for r = ResList
    ResAdj{r} = intersect(Reservoir(r).AdjacentRes,ResList);
end
colorRes = vertexcoloring(ResAdj,length(cmap(:,1)));

flowspac = 0; % 0.2; % spacing between flow lines

xlist = zeros(1,NumRes);
ylist = zeros(1,NumRes);

% Normalization of reservoir coordinates
x0 = Reservoir(ResList(1)).Centroid(1);
y0 = Reservoir(ResList(1)).Centroid(2);
x1 = Reservoir(Reservoir(ResList(1)).AdjacentRes(1)).Centroid(1);
y1 = Reservoir(Reservoir(ResList(1)).AdjacentRes(1)).Centroid(2);
dx0 = max([abs(x1 - x0) abs(y1 - y0)]);

hold on

if plotborders == 1
    % Plot the reservoirs
    xLinks = [];
    yLinks = [];
    i = 1;
    for r = ResList
        if ~isempty(Reservoir(r).BorderPoints)
            xr = coordscale(1).*(Reservoir(r).BorderPoints(1,:) - x0)./dx0;
            yr = coordscale(2).*(Reservoir(r).BorderPoints(2,:) - y0)./dx0;
            xLinks = [xLinks xr];
            yLinks = [yLinks yr];
            colori = cmap(colorRes(i),:);
            fill(xr,yr,colori,'EdgeColor','none')
            plot(xr,yr,'-','color',colori,'LineWidth',LW);
        end
        i = i + 1;
    end
    alpha(0.5);
end

for r = ResList
    xr = coordscale(1)*(Reservoir(r).Centroid(1) - x0)/dx0;
    yr = coordscale(2)*(Reservoir(r).Centroid(2) - y0)/dx0;
    xlist(r) = xr;
    ylist(r) = yr;
    
    % Plot flow exchanges
    if r < NumRes % to avoid flow line duplication
        for r2 = Reservoir(r).AdjacentRes
            if r2 > r && ismember(r2,ResList) % to avoid flow line duplication
                xj = coordscale(1)*(Reservoir(r2).Centroid(1) - x0)/dx0;
                yj = coordscale(2)*(Reservoir(r2).Centroid(2) - y0)/dx0;
                ang = atan2(yj-yr,xj-xr) + pi/2;
                dx = flowspac*cos(ang);
                dy = flowspac*sin(ang);
                xi1 = xr + dx;
                yi1 = yr + dy;
                xi2 = xr - dx;
                yi2 = yr - dy;
                xj1 = xj + dx;
                yj1 = yj + dy;
                xj2 = xj - dx;
                yj2 = yj - dy;
                % Effective flow from Ri to Rj
                plot([xi1 xj1],[yi1 yj1],'-','color','k','linewidth',LW)
                % Effective flow from Rj to Ri
                %plot([xi2 xj2],[yi2 yj2],'-','color','k','linewidth',LW)
            end
        end
    end
    
    % Plot reservoir disk
    th = 0:0.01:(2*pi);
    x = xr + ResRadius*cos(th);
    y = yr + ResRadius*sin(th);
    fill(x,y,'k','EdgeColor','none')
    
    text(xr,yr,['{\itR}_{' int2str(r) '}'],...
        'HorizontalAlignment','center','color',txtcolor,'FontName',fontname,'FontWeight','Bold','FontSize',FS)
end

% Plot size
xborder = 0.1; % increasing factor > 0 for the border spacing along x
yborder = 0.1; % increasing factor > 0 for the border spacing along x
if plotborders == 1
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
else
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
end

hold off
axis([xmin xmax ymin ymax])
daspect([1 1 1])
set(gca,'Position',[0 0 1 1],'FontName',fontname,'FontSize',FS)
set(gca,'visible','off')

end

