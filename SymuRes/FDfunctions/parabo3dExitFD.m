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

% Critical accumulation lines 
b1 = 2*bcc;
b2 = bbc + bcb;
b3 = vfc; % Critical acc. line for Global 3DMFD
      
% Global 3D MFD
ind_1 = find(b1*n(1,:) + b2*n(2,:) + b3 > 0); % Free flow side
ind_2 = find(b1*n(1,:) + b2*n(2,:) + b3 <= 0); % Congestion side

if ( ~isempty(ind_1) )
    Pt1 = parabo3dFD(n(:,ind_1),param);
    P(:,ind_1) = Pt1;
end
if ( ~isempty(ind_2) )
    Pt1 = parabo3dFD([max( -(b3 + b2*n(2,ind_2))/b1, 0 );n(2,ind_2)],param);
    P(:,ind_2) = Pt1;
end

end