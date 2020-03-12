function xres = findintersec2(xdata,ydata,yintersec)
% xres = findintersec2(xdata,ydata,yintersec)
% Find the intersection xres of the data (x,y) with a constant value y
% The data is supposed to define a monotonic function y = f(x)
%
% INPUTS
%---- xdata     : vector, abscissa of the data points (all different and
%                 sorted in ascending order)
%---- ydata     : vector, same size as xdata, ordinate of the data points
%---- yintersec : scalar, y value of the intersection
%
% OUTPUTS
%---- xres : scalar, result of the estimation f(xres) = y (linear interp.)

Ndata = length(xdata);

i = 1;
if ydata(1) < ydata(2) % if ydata is increasing
    while i <= Ndata && ydata(i) <= yintersec
        i = i + 1;
    end
else % if ydata is decreasing
    while i <= Ndata && ydata(i) >= yintersec
        i = i + 1;
    end
end
i = i - 1;

if i == 0
    xres = -Inf; % the value of yintersec is too low, no intersection found
elseif i == Ndata
    xres = +Inf; % the value of yintersec is too high, no intersection found
else
    % the abscissa of the intersection is estimated by linear interpolation
    xres = xdata(i) + (xdata(i+1) - xdata(i))/(ydata(i+1) - ydata(i))*(yintersec - ydata(i));
end

end