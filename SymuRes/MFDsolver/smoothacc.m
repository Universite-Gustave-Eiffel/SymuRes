function nsmooth = smoothacc(n0,tin,tout,Nsmooth)
% nsmooth = smoothacc(n0,tin,tout,Nsmooth)
% Smooth the variation of accumulation in calculating the exact non integer
% value of n thanks to the count curves
%
% INPUTS
%---- n0      : integer, actual integer value of accumulation
%---- tin     : vector, list of entry times
%---- tout    : vector, list of exit times
%---- Nsmooth : integer, number of previous times used to approximate the
%               local evolution of the count curves
%
% OUTPUTS
%---- nsmooth : scalar, approximated (smoothed) value of accumulation

if length(tout) > Nsmooth + 1
    if tin(end) > tout(end) % the last event was an entry
        t1 = tin(end);
        t2 = tout(end);
        N1 = length(tin) - 1;
        N2 = length(tout) - 1;
        N_out = 1:length(tout);
        % estimation of the local qout by linear regression
        q1 = linregr2(tout((end-Nsmooth):end),N_out((end-Nsmooth):end));
        nsmooth = N1 - N2 - (t1 - t2)*q1;
    else % the last event was an exit
        t1 = tout(end);
        t2 = tin(end);
        N1 = length(tout) - 1;
        N2 = length(tin) - 1;
        N_in = 1:length(tin);
        % estimation of the local qin by linear regression
        q1 = linregr2(tin((end-Nsmooth):end),N_in((end-Nsmooth):end));
        nsmooth = N2 - N1 + (t1 - t2)*q1;
    end
else
    nsmooth = n0;
end

end