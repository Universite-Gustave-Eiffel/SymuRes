function singlearrow(x,y,arrowLfactor,arrowLtype,arrowWfactor,LColor,LW)
% singlearrow(x,y,arrowL,arrowW)
% Plot a single arrow
%
% INPUTS
%---- x, y         : vectors, abscissa and ordinate of the extremities of the
%                    arrow [x(1) x(2)], [y(1) y(2)]
%---- arrowLfactor : scalar > 0, arrow length
%---- arrowLtype   : string, arrow length value, 'absolute' or 'relative'
%                    (i.e. length as the ratio over the line length)
%---- arrowWfactor : scalar > 0, arrow width ratio over the arrow length
%---- LColor       : color, string or RGB vector
%---- LW           : scalar, arrow line width

dasp = daspect; % aspect ratio of the current axes
xratio = dasp(1);
yratio = dasp(2);
xscale = xratio/max(xratio,yratio);
yscale = yratio/max(xratio,yratio);
% xscale = 1;
% yscale = yratio/xratio;

hold on

% Plot the line
plot(x,y,'-','color',LColor,'linewidth',LW)

% Plot arrow
i = 1;
alpha = atan2(1/yscale*(y(i+1) - y(i)),1/xscale*(x(i+1) - x(i))); % line angle
if strcmp(arrowLtype,'relative')
    L = sqrt(((x(i+1) - x(i))/xscale)^2 + ((y(i+1) - y(i))/yscale)^2); % line length
    arrowL = arrowLfactor*L;
else
    arrowL = arrowLfactor;
end
arrowW = arrowWfactor*arrowL;
x1 = x(i+1) - xscale*(arrowL*cos(alpha) + arrowW/2*sin(alpha));
y1 = y(i+1) - yscale*(arrowL*sin(alpha) - arrowW/2*cos(alpha));
x2 = x(i+1) - xscale*(arrowL*cos(alpha) - arrowW/2*sin(alpha));
y2 = y(i+1) - yscale*(arrowL*sin(alpha) + arrowW/2*cos(alpha));
plot([x(i+1) x1],[y(i+1) y1],'-',[x(i+1) x2],[y(i+1) y2],'-','color',LColor,'linewidth',LW)

% th = 0:0.01:2*pi;
% plot(x(i)+xscale*arrowL*cos(th),y(i)+yscale*arrowL*sin(th),'-','color',LColor,'linewidth',LW)

end