function P = parabo3dExitFD(n,param)
% P = parabo3dExitFD(n,param)
% Compute the 3D parabolic exit function (from the parabolic FD model)
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
%---- P : matrix [P_c; P_b], car and bus exit demand productions [veh.m/s]

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

% Car 3D MFD
% P_c = n_c * (vfc + bcc * n_c + bbc * n_b)
% dP_c/dn_c = vfc + 2 * bcc * n_c + bbc * n_b
% To maximum P_c, dP_c/dn_c = 0 => vfc + 2 * bcc * n_c + bbc * n_b = 0
% For maximum P_c for a given n_c, n_b should be (vfc - bbc * n_b) / (2 * bcc)
ind_1 = find(2*bcc*n(1,:) + bbc*n(2,:) + vfc <= 0); % Congestion side
ind_2 = find(2*bcc*n(1,:) + bbc*n(2,:) + vfc > 0); % Free flow side

if ( ~isempty(ind_1) )
    Pt = parabo3dFD([max( -(vfc + bbc*n(2,ind_1))/(2*bcc), 0 ); n(2,ind_1)],param);
    P(1,ind_1) = Pt(1);
end
if ( ~isempty(ind_2) )
    Pt = parabo3dFD(n(:,ind_2),param);
    P(1,ind_2) =  Pt(1);
end

% Bus 3D MFD
% P_b = n_b * (vfb + bcb * n_c + bbb * n_b)
% dP_b/dn_b = vfb + bcb * n_c + 2 * bbb * n_b
% To maximum P_b, dP_b/dn_b = 0 => vfb + bcb * n_c + 2 * bbb * n_b = 0
% For maximum P_b for a given n_b, n_c should be (vfb - bcb * n_b) / (2 * bbb)
ind_1 = find(bcb*n(1,:) + 2*bbb*n(2,:) + vfb <= 0); % Congestion side
ind_2 = find(bcb*n(1,:) + 2*bbb*n(2,:) + vfb > 0); % Free flow side

if ( ~isempty(ind_1) )
    Pt = parabo3dFD([n(1,ind_1);max( -(vfb + bcb*n(1,ind_1))/(2*bbb), 0 )],param);
    P(2,ind_1) = Pt(2);
end
if ( ~isempty(ind_2) )
    Pt = parabo3dFD(n(:,ind_2),param);
    P(2,ind_2) =  Pt(2);
end

end