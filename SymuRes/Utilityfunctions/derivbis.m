function sderiv = derivbis(t,signal)

Nsignal = length(signal);
sderiv = zeros(1,Nsignal);

for i = 1:(Nsignal-1)
    sderiv(i) = (signal(i+1) - signal(i))/(t(i+1) - t(i));
end
sderiv(Nsignal) = sderiv(Nsignal-1);

end