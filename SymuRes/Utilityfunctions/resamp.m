function xsample = resamp(tsample,tdata,xdata)
% xsample = resamp(tsample,tdata,xdata)
% Return data signal points resampled according to a new time
% discretization, using linear interpolation.
%
% INPUTS
%---- tsample : row vector, new time sample
%---- tdata   : row vector, old time discretization, may be non regular
%---- xdata   : row vector, same size as tdata, signal points at each time
%
% OUTPUTS
%---- xsample : row vector, same size as tsample, signal points resampled
%               at the new time points (linear interpolation)

Ns = length(tsample);
Nt = length(tdata);

xsample = zeros(1,Ns);

i = 1;
for j = 1:Ns
    while i <= Nt && tdata(i) <= tsample(j)
        i = i + 1;
    end
    i = max([i-1 1]);
    
    if i < Nt
        xsample(j) = xdata(i) + (tsample(j) - tdata(i))*(xdata(i+1) - xdata(i))/(tdata(i+1) - tdata(i));
    else
        xsample(j) = xdata(i);
    end
end

end 