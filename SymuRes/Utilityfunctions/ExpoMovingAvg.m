function savg = ExpoMovingAvg(signal,smoothdegree)
% savg = ExpoMovingAvg(signal,smoothdegree)
% Compute the Exponential Moving Average of a signal
% 
% INPUTS
%---- signal       : row vector, signal values
%---- smoothdegree : scalar in [0,1], small value gives highly smoothed results
%
% OUTPUTS
%---- savg : row vector (same size as signal), resulting "smoothed" signal

Nsignal = length(signal);
savg = zeros(1,Nsignal);

savg(1) = signal(1);
for i = 2:Nsignal
    savg(i) = savg(i-1) + smoothdegree*(signal(i) - savg(i-1));
end

end