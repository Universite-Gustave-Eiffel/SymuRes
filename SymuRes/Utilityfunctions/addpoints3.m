function [x, y, t] = addpoints3(x0,y0,t0,dtmin)
% [x, y, t] = addpoints3(x0,y0,t0,dtmin)
% Add points regularly spaced between each point of the data (x,y,t) so
% that there is at least a time diff of tmin between two consecutive points
%
% INPUTS
%---- x0,y0,t0 : row vectors, original data of abscissa, ordinate and time
%---- dtmin    : scalar, minimum time diff required between two points
%
% OUTPUTS
%---- x,y,t : row vectors, data artificially completed with linear interpolations

Ndata = length(x0);
Ndata2 = 10*Ndata;

x = zeros(1,Ndata2);
y = zeros(1,Ndata2);
t = zeros(1,Ndata2);

j0 = 1;
for i = 1:(Ndata-1)
    dt = t0(i+1) - t0(i); % distance between pt i and i+1
    Npts = floor(dt/dtmin); % number of intermediate points to add
    
    x(j0) = x0(i);
    y(j0) = y0(i);
    t(j0) = t0(i);
    if Npts > 0
        for j = 1:Npts
            x(j0+j) = x0(i) + j/(Npts+1)*(x0(i+1) - x0(i));
            y(j0+j) = y0(i) + j/(Npts+1)*(y0(i+1) - y0(i));
            t(j0+j) = t0(i) + j/(Npts+1)*(t0(i+1) - t0(i));
        end
    end
    j0 = j0 + 1 + Npts;
end
x(j0) = x0(Ndata);
y(j0) = y0(Ndata);
t(j0) = t0(Ndata);
x = x(1:j0);
y = y(1:j0);
t = t(1:j0);

end