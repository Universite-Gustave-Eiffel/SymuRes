function savg = SimpleMovingAvg(signal,period)
% savg = SimpleMovingAvg(signal,period)
% Compute the Simple Centered Moving Average of a signal with a given
% period
% 
% INPUTS
%---- signal : row vector, signal values
%---- period : integer, number of consecutive signal values on which the arithmetic mean is calculated
%
% OUTPUTS
%---- savg : row vector (same size as signal), resulting "smoothed" signal

Nsignal = length(signal);
savg = zeros(1,Nsignal);

if floor(period/2) == period/2
    period2 = period/2;
    for i = (period2 + 1):(Nsignal - period2 + 1)
        savg(i) = sum(signal((i - period2):(i + period2 - 1)))/period;
    end
    for i = 1:(period2 + 1)
        savg(i) = savg(period2 + 1);
    end
    for i = (Nsignal - period2 + 1):Nsignal
        savg(i) = savg(Nsignal - period2 + 1);
    end
else
    period2 = floor(period/2);
    for i = (period2 + 1):(Nsignal - period2)
        savg(i) = sum(signal((i - period2):(i + period2)))/period;
    end
    for i = 1:(period2 + 1)
        savg(i) = savg(period2 + 1);
    end
    for i = (Nsignal - period2):Nsignal
        savg(i) = savg(Nsignal - period2);
    end
end

end