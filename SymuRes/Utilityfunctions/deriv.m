function sderiv = deriv(t,signal)
% sderiv = deriv(t,signal)
% Return the first derivative of a signal, respect to time
%
% INPUTS
%---- t      : vector, times
%---- signal : vector, same size as t, signal values
%
% OUTPUTS
%---- sderiv : vector, same size as t, derivative of the signal

Nsignal = length(signal);
sderiv = zeros(1,Nsignal);

for i = 2:(Nsignal-1)
    sa = (signal(i) - signal(i-1))/(t(i) - t(i-1));
    sb = (signal(i+1) - signal(i))/(t(i+1) - t(i));
    sderiv(i) = ((t(i+1) - t(i))*sa + (t(i) - t(i-1))*sb)/(t(i+1) - t(i-1));
end

i = 2;
sa = (signal(i) - signal(i-1))/(t(i) - t(i-1));
sderiv(1) = 2*sa - sderiv(i);

i = Nsignal-1;
sb = (signal(i+1) - signal(i))/(t(i+1) - t(i));
sderiv(Nsignal) = 2*sb - sderiv(i);

end