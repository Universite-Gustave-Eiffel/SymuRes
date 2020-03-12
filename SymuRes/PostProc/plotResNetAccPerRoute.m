function plotResNetAccPerRoute(t,Link,Reservoir,ResList,Route,RoutesList,SimulTime,opts)
% plotResNetAccPerRoute(t,Link,Reservoir,ResList,Route,RoutesList,CurrentTime,opts)
% Plot the state of reservoirs at time t (accumulation), with links and/or shape borders
% Plot accumulation ratio of each route in the reservoirs
%
% INPUTS
%---- t           : scalar, time [s]
%---- Link        : Link structure, put [] if no links to plot
%---- Reservoir   : Reservoir structure
%---- ResList     : vector, reservoir IDs
%---- Route       : Route structure
%---- RoutesList  : vector, route IDs
%---- SimulTime   : vector, simulation time [s]
%---- opts        : options, structure with fields 'fontname', 'fontsize',
%                   'linewidth', 'colormap', 'textcolor', 'title', 'showleg', 'legloc'

NbL = length(Link);
NbR = length(ResList);
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
    FS = 14; % default
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
if isfield(opts,'textcolor')
    txtcolor = opts.textcolor;
else
    txtcolor = [204 0 0]/255; % default
end
if isfield(opts,'title')
    tlabel = opts.title;
else
    tlabel = {['\bf{\itt} = ' num2str(t) ' s']; '\rmAccumulation [veh]'}; % default
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
colorBorder = [0.1 0.1 0];
hf = zeros(1,NbRoutes);
strleg = cellstr(int2str(zeros(NbRoutes,1)));
ileg = 1;
LegList = [];
for r = ResList
    xr = Reservoir(r).Centroid(1);
    yr = Reservoir(r).Centroid(2);
    xlist = Reservoir(r).BorderPoints(1,:);
    ylist = Reservoir(r).BorderPoints(2,:);
    
    dist = sqrt((xr - xlist).^2 + (yr - ylist).^2);
    maxdist = 1.01*max(dist);
    
    angi = 0;
    xstart = xr + maxdist*cos(angi);
    ystart = yr + maxdist*sin(angi);
    i = 1;
    isinter = 0;
    while isinter == 0 && i < length(xlist)
        [isinter, xb, yb] = lineintersect([xr yr],[xstart ystart],[xlist(i) ylist(i)],[xlist(i+1) ylist(i+1)]);
        i = i + 1;
    end
    indstart = i - 1;
    xborderstart = xb;
    yborderstart = yb;
    indstart0 = i - 1;
    xborderstart0 = xb;
    yborderstart0 = yb;
    
    k_r = 1;
    for iroute = RoutesList
        i_r = find(Reservoir(r).RoutesID == iroute);
        if ~isempty(i_r)
            accratio = Reservoir(r).AccPerRoute(i_r,timeID)/Reservoir(r).MaxAcc;
            angi = angi + accratio*2*pi;
            xend = xr + maxdist*cos(angi);
            yend = yr + maxdist*sin(angi);
            indend = indstart;
            for i = [indstart:length(xlist) 1:(indstart-1)]
                if i == length(xlist)
                    [isinter, xb, yb] = lineintersect([xr yr],[xend yend],[xlist(i) ylist(i)],[xlist(1) ylist(1)]);
                else
                    [isinter, xb, yb] = lineintersect([xr yr],[xend yend],[xlist(i) ylist(i)],[xlist(i+1) ylist(i+1)]);
                end
                if isinter == 1
                    indend = i;
                    xborderend = xb;
                    yborderend = yb;
                end
            end
            if indstart <= indend
                xlist0 = xlist(indstart:indend);
                ylist0 = ylist(indstart:indend);
            else
                xlist0 = [xlist(indstart:end) xlist(1:indend)];
                ylist0 = [ylist(indstart:end) ylist(1:indend)];
            end
            if ~ismember(iroute,LegList) % add to the legend
                hf(ileg) = fill([xr xborderstart xlist0 xborderend],[yr yborderstart ylist0 yborderend],...
                    cmap(k_r,:),'EdgeColor','none');
                strleg{ileg} = ['[' int2str(Route(iroute).ResPath) ']'];
                LegList = [LegList iroute];
                ileg = ileg + 1;
            else
                fill([xr xborderstart xlist0 xborderend],[yr yborderstart ylist0 yborderend],...
                    cmap(k_r,:),'EdgeColor','none');
            end
            indstart = indend;
            xborderstart = xborderend;
            yborderstart = yborderend;
        end
        k_r = k_r + 1;
    end
    
    indend = indstart0;
    xborderend = xborderstart0;
    yborderend = yborderstart0;
    if indstart < indend
        xlist0 = xlist(indstart:indend);
        ylist0 = ylist(indstart:indend);
    else
        xlist0 = [xlist(indstart:end) xlist(1:indend)];
        ylist0 = [ylist(indstart:end) ylist(1:indend)];
    end
    plot(Reservoir(r).BorderPoints(1,:),Reservoir(r).BorderPoints(2,:),'-','color',colorBorder,'LineWidth',LW);
end
for r = ResList
    xr = Reservoir(r).Centroid(1);
    yr = Reservoir(r).Centroid(2);
    text(xr,yr,{['{\itR}_{' int2str(r) '}']; int2str(round(Reservoir(r).Acc(timeID)))},...
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

text((xmin+xmax)/2,ymax-yborder*dy/2,tlabel,...
    'HorizontalAlignment','center','FontName',fontname,'FontSize',FS)
hold off
xlabel('\itx \rm[m]','FontName',fontname,'FontSize',FS)
ylabel('\ity \rm[m]','FontName',fontname,'FontSize',FS)
set(gca,'Position',[0 0 1 1],'FontName',fontname,'FontSize',FS)
set(gca,'visible','off')
if showleg == 1
    hleg = legend(hf,strleg);
    set(hleg,'Location',legloc,'FontName',fontname,'FontSize',0.8*FS)
end

end