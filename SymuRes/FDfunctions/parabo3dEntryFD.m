function P = parabo3dEntryFD(n,param)
% P = parabo3dEntryFD(n,param)
% Compute the 3D parabolic entry function (from the parabolic FD model)
%
% Feb 2020 - Mahendra Paipuri
%
% INPUTS
%---- n     : matrix [n_c; n_b], car and bus accumulations [veh]
%---- param : row vector, parameters of the FD function: [vfc vfb bcc bbc bcb bbb nc0 nb0]
%             -vfc: scalar, free flow speed of cars [m/s]
%             -vfb: scalar, free flow speed of buses [m/s]
%             -bcc: scalar, marginal effect of cars on free flow speed of cars [m/s/veh]
%             -bbc: scalar, marginal effect of buses on free flow speed of cars [m/s/veh]
%             -bcb: scalar, marginal effect of cars on free flow speed of buses [m/s/veh]
%             -bbb: scalar, marginal effect of buses on free flow speed of buses [m/s/veh]
%             -nc0: scalar, normalization accumulation for cars, default to 1 [veh]
%             -nb0: scalar, normalization accumulation for buses, default to 1 [veh]
%
% OUTPUTS
%---- P : matrix [P_c; P_b], car and bus entry supply productions [veh.m/s]

if ( size(n,1) ) == 1
    n = [n(1,:); zeros(size(n(1,:)))];    
end

nlen = size(n,2);
P = zeros(2,nlen);

vfc = param(1)/param(7);
vfb = param(2)/param(8);
bcc = param(3)/param(7)^2;
bbc = param(4)/param(7)^2;
bcb = param(5)/param(8)^2;
bbb = param(6)/param(8)^2;

     
% Global 3D MFD
% P = n_c * (vfc + bcc * n_c + bbc * n_b) + n_b * (vfb + bcb * n_c + bbb * n_b)
% dP/dn_c = vfc + 2 * bcc * n_c + (bbc + bcb) * n_b
% To maximum P, dP/dn_c = 0 => vfc + 2 * bcc * n_c + (bbc + bcb) * n_b = 0
% For maximum P for a given n_c, n_b should be (vfc - ((bbc + bcb)) * n_b) / (2 * bcc)
cong_ind = find(2*bcc*n(1,:) + (bbc + bcb)*n(2,:) + vfc <= 0); % Congestion side
ff_ind = find(2*bcc*n(1,:) + (bbc + bcb)*n(2,:) + vfc > 0); % Free flow side

if ( ~isempty(cong_ind) )
    Pt = parabo3dFD(n(:,cong_ind),param);
    P(:,cong_ind) = Pt;
end
if ( ~isempty(ff_ind) )
    Pt = parabo3dFD([max( -(vfc + (bbc + bcb)*n(2,ff_ind))/(2*bcc), 0 ); n(2,ff_ind)],param);
    P(:,ff_ind) =  Pt;
end

% Bus 3D MFD
% P_b = n_b * (vfb + bcb * n_c + bbb * n_b)
% dP_b/dn_b = vfb + bcb * n_c + 2 * bbb * n_b
% To maximum P_b, dP_b/dn_b = 0 => vfb + bcb * n_c + 2 * bbb * n_b = 0
% For maximum P_b for a given n_b, n_c should be (vfb - bcb * n_b) / (2 * bbb)
% ind_1 = find(bcb*n(1,:) + 2*bbb*n(2,:) + vfb <= 0); % Congestion side
% ind_2 = find(bcb*n(1,:) + 2*bbb*n(2,:) + vfb > 0); % Free flow side
% 
% if ( ~isempty(ind_1) )
%     Pt = parabo3dFD(n(:,ind_1),param);
%     P(2,ind_1) = Pt(2);
% end
% if ( ~isempty(ind_2) )
%     Pt = parabo3dFD([n(1,ind_2);max( -(vfb + bcb*n(1,ind_2))/(2*bbb), 0 )],param);
%     P(2,ind_2) =  Pt(2);
% end

end