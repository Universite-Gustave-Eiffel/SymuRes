function [x, y] = addpoints(x0,y0,Npts)
% [x, y] = addpoints(x0,y0,Npts)
% Add Npts points regularly spaced between each point of the data (x,y)
%
% INPUTS
%---- x0,y0 : row vectors, original data
%---- Npts  : integer, number of points to add between each original point
%
% OUTPUTS
%---- x,y : row vectors, data artificially completed with linear interpolations

Ndata = length(x0);
Ndata2 = (Npts+1)*(Ndata-1) + 1;

x = zeros(1,Ndata2);
y = zeros(1,Ndata2);

for i = 1:(Ndata-1)
    x(1+(Npts+1)*(i-1)) = x0(i);
    y(1+(Npts+1)*(i-1)) = y0(i);
    for j = 1:Npts
        x(j+1+(Npts+1)*(i-1)) = x0(i) + 1/(Npts+1)*(x0(i+1) - x0(i));
        y(j+1+(Npts+1)*(i-1)) = y0(i) + 1/(Npts+1)*(y0(i+1) - y0(i));
    end
end
x(Ndata2) = x0(Ndata);
y(Ndata2) = y0(Ndata);

end