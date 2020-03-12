function [a, b, r, da, xregr, slopemargin, datamargin] = linregr(X,Y,alpha)
% [a, b, r, da, xregr, slopemargin, datamargin] = linregr(X,Y,alpha)
% Compute simple linear regression of the data Y explained by data X (each
% point in X has a unique corresponding point in Y): Y = a.X + b
%
% INPUTS
%---- X,Y   : vectors or matrix, contains the data to analyze 
%---- alpha : scalar in [0,1], risk parameter of the Student's law
%
% OUTPUT
%---- a,b         : scalars, estimation of the linear regression parameters
%---- r           : scalar, coefficient of correlation
%---- da          : scalar, confidence margin on the estimation of the slope a
%---- xregr       : row vector, sample data (equally space between xmin and
%                   xmax) to plot the regression line
%---- slopemargin : row vector (same size as xregr), confidence margin on
%                   estimated points of the regression
%---- datamargin  : row vector (same size as xregr), confidence margin on
%                   data

n = size(Y,1)*size(Y,2);
ddl = n - 2; % degree of freedom
xm = mean(mean(X));
xse2 = sum(sum((X - xm).^2));
ym = mean(mean(Y));
yse2 = sum(sum((Y - ym).^2));
xyse = sum(sum((X - xm).*(Y - ym)));

a = xyse/xse2; % slope estimation
b = ym - a*xm; % y-intercept

r = a*sqrt(xse2/yse2); % coefficient of correlation

stddev = sqrt((1 - r^2)/ddl*yse2); % standard deviation for the error on data
ta2 = abs(tstudent(alpha/2,ddl)); % variable of the Student's law for risk alpha/2
da = ta2*stddev/xse2; % confidence margin on the estimation of the slope a

xregr = linspace(min(min(X)),max(max(X)),100);
slopemargin = ta2*stddev*sqrt(1/n + (xregr - xm).^2./xse2);
datamargin = ta2*stddev*sqrt(1 + 1/n + (xregr - xm).^2./xse2);

end