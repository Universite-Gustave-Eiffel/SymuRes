function singlearrowtext(x,y,strtxt,arrowL,arrowW,LColor,LW,fontname,FS,horzalign,vertalign)
% singlearrowtext(x,y,strtxt,arrowL,arrowW,LColor,LW,fontname,FS,horzalign,vertalign)
% Plot a single arrow with a text at the second extremity
%
% INPUTS
%---- x, y      : vectors, abscissa and ordinate of the extremities of the
%                 arrow [x(1) x(2)], [y(1) y(2)]
%---- strtxt    : string, text to print at the second extremity of the arrow
%---- arrowL    : scalar, arrow length
%---- arrowW    : scalar, arrow width
%---- LColor    : color, string or RGB vector
%---- LW        : scalar, arrow line width
%---- fontname  : string, font name
%---- FS        : integer, font size
%---- horzalign : string, text horizontal alignment (see text properties)
%---- vertalign : string, text vertical alignment (see text properties)

n = length(x);

dasp = daspect; % aspect ratio of the current axes
xratio = dasp(1);
yratio = dasp(2);
xscale = xratio/max(xratio,yratio);
yscale = yratio/max(xratio,yratio);

hold on

% Plot the line
plot(x,y,'-','color',LColor,'linewidth',LW)

% Plot the first arrow
i = 1;
alpha = atan2(1/yscale*(y(i) - y(i+1)),1/xscale*(x(i) - x(i+1)));
x1 = x(i) - xscale*(arrowL*cos(alpha) + arrowW/2*sin(alpha));
y1 = y(i) - yscale*(arrowL*sin(alpha) - arrowW/2*cos(alpha));
x2 = x(i) - xscale*(arrowL*cos(alpha) - arrowW/2*sin(alpha));
y2 = y(i) - yscale*(arrowL*sin(alpha) + arrowW/2*cos(alpha));
plot([x(i) x1],[y(i) y1],'-',[x(i) x2],[y(i) y2],'-','color',LColor,'linewidth',LW)

% Plot text
text(x(2),y(2),strtxt,'HorizontalAlignment',horzalign,'VerticalAlignment',vertalign,'color',LColor,'FontName',fontname,'FontSize',FS)

end