function plotResNetSpeed(t,Link,Reservoir,ResList,SimulTime,SpeedRange,colorscheme,opts)
% plotResNetSpeed(t,Link,Reservoir,ResList,CurrentTime,SpeedRange,colorscheme,opts)
% Plot the state of reservoirs at time t (mean speed), with links and/or shape borders
%
% INPUTS
%---- t           : scalar, time [s]
%---- Link        : Link structure, put [] if no links to plot
%---- Reservoir   : Reservoir structure
%---- ResList     : vector, reservoir IDs
%---- SimulTime   : vector, simulation time [s]
%---- SpeedRange  : vector [Vmin Vmax], speed range [m/s] to define the colormap
%---- colorscheme : string, color rendering 'jet' or 'trafficolor'
%---- opts        : options, structure with fields 'fontname', 'fontsize',
%                   'linewidth', 'title', 'legloc'

% timeID = floor(t/TimeStep) + 1; % index of the current time
timeID = findindex(SimulTime,t); % index of the current time

% Choice of a colormap
nbColor = 800;
if strcmp(colorscheme,'jet')
    cmap = colormap(jet(nbColor)); % jet colormap
    txtcolor = [0.9 0.9 1];
elseif strcmp(colorscheme,'trafficolor')
    freeflowColor = [0 1 0]; % green
    congestionColor = [1 0 0]; % red
    cmap = buildmap(freeflowColor,congestionColor,nbColor); % trafficolor colormap (green-red)
    %txtcolor = [0.9 0.9 1];
    txtcolor = [0.1 0.1 0];
end
% cmap(1,:) = [0 0 0]; % black color for null density
cmap = reversematrix(cmap,1);
colormap(cmap);

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
if isfield(opts,'linewidth')
    LW = opts.linewidth;
else
    LW = 2; % default
end
if isfield(opts,'title')
    tlabel = opts.title;
else
    tlabel = ['\itt \rm\bf= ' num2str(t) ' s']; % default
end
if isfield(opts,'legloc')
    legloc = opts.legloc;
else
    legloc = 'EastOutside'; % default
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
        xLinks = [xLinks Reservoir(r).BorderPoints(1,:)];
        yLinks = [yLinks Reservoir(r).BorderPoints(2,:)];
    end
end

% Plot the reservoirs
vmax = 0;
for r = ResList
    if Reservoir(r).FreeflowSpeed > vmax
        vmax = Reservoir(r).FreeflowSpeed;
    end
end
for r = ResList
    %accratio = Reservoir(r).Acc(timeID)/Reservoir(r).MaxAcc;
    %speedratio = Reservoir(r).MeanSpeed(timeID)/vmax;
    speedratio = (Reservoir(r).MeanSpeed(timeID) - SpeedRange(1))/(SpeedRange(2) - SpeedRange(1));
    indcolor = min([max([floor(speedratio*nbColor) 1]) nbColor]);
    colori = cmap(indcolor,:);
    
    fill(Reservoir(r).BorderPoints(1,:),Reservoir(r).BorderPoints(2,:),colori,'EdgeColor','none')
    plot(Reservoir(r).BorderPoints(1,:),Reservoir(r).BorderPoints(2,:),'-','color',colori,'LineWidth',LW);
end
for r = ResList
    xr = Reservoir(r).Centroid(1);
    yr = Reservoir(r).Centroid(2);
    text(xr,yr,{['{\itR}_{' int2str(r) '}']; [int2str(round(Reservoir(r).MeanSpeed(timeID)*3.6)) ' km/h']},...
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

hcb = colorbar(legloc);
hcb.Label.String = 'Mean speed [km/h]';
hcb.Label.FontSize = FS;
hcb.Label.FontName = fontname;
caxis(3.6*SpeedRange);
axis([xmin xmax ymin ymax])
daspect([1 1 1])
alpha(0.5);

text((xmin+xmax)/2,ymax-yborder*dy/4,tlabel,...
    'HorizontalAlignment','center','FontName',fontname,'FontSize',FS)
xlabel('\itx \rm[m]','FontName',fontname,'FontSize',FS)
ylabel('\ity \rm[m]','FontName',fontname,'FontSize',FS)
set(gca,'Position',[0 0.02 1 0.96],'FontName',fontname,'FontSize',FS)
set(gca,'visible','off')
% set(gcf,'Position',[10 10 1000 800])
hold off

end