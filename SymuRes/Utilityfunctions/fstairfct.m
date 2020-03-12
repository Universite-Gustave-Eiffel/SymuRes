function fval = fstairfct(t,f,Nstairs)
% fval = fstairfct(t,f,Nstairs)
% Return a function approximated by stairs
%
% INPUTS
%---- t       : scalar or vector, time [s]
%---- f       : function handle
%---- Nstairs : scalar, number of stairs
%
% OUTPUTS
%--- fval : vector, same size as t, stair function

Nt = length(t);
Npts = 1000; % number of points to calculate the mean value for a stair

timestair = linspace(t(1),t(Nt),Nstairs+1); % time of each stair edges

fval = zeros(1,Nt);

for i = 1:(Nstairs-1)
    stairt = linspace(timestair(i),timestair(i+1),Npts);
    stairvalue = mean(f(stairt));
    fval = fval + (timestair(i) <= t).*(t < timestair(i+1)).*stairvalue;
end

i = Nstairs;
stairt = linspace(timestair(i),timestair(i+1),Npts);
stairvalue = mean(f(stairt));
fval = fval + (timestair(i) <= t).*(t <= timestair(i+1)).*stairvalue;

end