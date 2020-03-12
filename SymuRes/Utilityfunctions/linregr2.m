function [a, b] = linregr2(x,y)
% [a, b] = linregr2(x,y)
% Compute simple linear regression of the data y explained by data x (each
% point in x has a unique corresponding point in y): y = a.x + b
%
% INPUTS
%---- x,y : vectors of same size, contains the data to analyze 
%
% OUTPUT
%---- a,b : scalars, estimation of the linear regression parameters

xm = mean(x);
xse2 = sum((x - xm).^2);
ym = mean(y);
xyse = sum((x - xm).*(y - ym));

a = xyse/xse2; % slope estimation
b = ym - a*xm; % y-intercept

end