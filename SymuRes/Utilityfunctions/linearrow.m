function linearrow(x,y,ngap,arrowL,arrowW,LColor,LW)
% linearrow(x,y,ngap,arrowL,arrowW,LColor,LW)
% Plot arrows along a 2D curve to indicate the direction
%
% INPUTS
%---- x, y   : vectors, abscissa and ordinate of the points of the curve
%---- ngap   : scalar, plot an arrow each ngap points
%---- arrowL : scalar, arrow length
%---- arrowW : scalar, arrow width
%---- LColor : color, string or RGB vector
%---- LW     : scalar, arrow line width

n = length(x);

dasp = daspect; % aspect ratio of the current axes
xratio = dasp(1);
yratio = dasp(2);
xscale = xratio/max(xratio,yratio);
yscale = yratio/max(xratio,yratio);

hold on
for i=ngap:ngap:n
    alpha = atan2(1/yscale*(y(i) - y(i-1)),1/xscale*(x(i) - x(i-1)));
    x1 = x(i) - xscale*(arrowL*cos(alpha) + arrowW/2*sin(alpha));
    y1 = y(i) - yscale*(arrowL*sin(alpha) - arrowW/2*cos(alpha));
    x2 = x(i) - xscale*(arrowL*cos(alpha) - arrowW/2*sin(alpha));
    y2 = y(i) - yscale*(arrowL*sin(alpha) + arrowW/2*cos(alpha));
    
    plot([x(i) x1],[y(i) y1],'-',[x(i) x2],[y(i) y2],'-','color',LColor,'linewidth',LW)
end

end