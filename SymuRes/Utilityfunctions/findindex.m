function ind0 = findindex(xdata,x0)
% ind = findindex(xdata,x0)
% Find the index of the point in xdata closest to x0
%
% INPUTS
%---- xdata : vector, length > 2, values must be sorted in ascending order
%---- x0    : scalar or vector, point for which we want to find the closest value
%
% OUTPUTS
%---- ind0 : vector of integers, same size as x0, index of the closest value in xdata

Ndata = length(xdata);
Nx0 = length(x0);

ind0 = zeros(1,Nx0);

for j = 1:Nx0
    xj = x0(j);
    for i = 2:(Ndata-1)
        x1 = (xdata(i) + xdata(i-1))/2;
        x2 = (xdata(i+1) + xdata(i))/2;
        if x1 <= xj && xj < x2
            ind = i;
        end
    end
    
    i = 1;
    x2 = (xdata(i+1) + xdata(i))/2;
    if xj < x2
        ind = i;
    end
    
    i = Ndata;
    x1 = (xdata(i) + xdata(i-1))/2;
    if x1 <= xj
        ind = i;
    end
    ind0(j) = ind;
end

end