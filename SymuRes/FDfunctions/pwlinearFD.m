function q = pwlinearFD(k,param)
% q = pwlinearFD(k,param)
% Compute the piecewise linear Fundamental Diagram function
%
% INPUTS
%---- k     : scalar or vector, vehicle density [veh/m], range of k must be between 0 and kj
%---- param : Npts-by-2 matrix, points to define the piecewise linear FD
%             column 1 contains the densities in increasing order [veh/m]
%             column 2 contains the flows corresponding to each density [veh/s]
%
% OUTPUTS
%---- q : scalar or vector (same size as k), vehicle flow [veh/s]

kj = param(end,1); % jam density (max. density) [veh/m]
pts = param(1:(end-1),:);

Npts = size(pts,1); % number of points for the definition of the FD
kc = pts(:,1);
qc = pts(:,2);
q = zeros(1,length(k));

if kc == sort(kc) % test if kc is sorted by increasing order
    kc0 = [0 kc' kj];
    qc0 = [0 qc' 0];
    dqc0 = derivbis(kc0,qc0);
    if dqc0 == sort(dqc0,'descend') % test if qc really gives a concave FD
        if Npts == 1
            kc1 = kc(1);
            qc1 = qc(1);
            q = (k <= kc1).*(qc1/kc1.*k) + (kc1 < k).*(k < kj).*(qc1 - qc1/(kj - kc1).*(k - kc1));
        else
            for i = 1:(Npts-1)
                kc1 = kc(i);
                kc2 = kc(i+1);
                qc1 = qc(i);
                qc2 = qc(i+1);
                q = q + (kc1 < k).*(k <= kc2).*(qc1 + (qc2 - qc1)/(kc2 - kc1).*(k - kc1));
            end
            kc1 = kc(1);
            qc1 = qc(1);
            q = q + (k <= kc1).*(qc1/kc1.*k);
            kc2 = kc(Npts);
            qc2 = qc(Npts);
            q = q + (kc2 < k).*(k < kj).*(qc2 - qc2/(kj - kc2).*(k - kc2));
        end
    end
end

end