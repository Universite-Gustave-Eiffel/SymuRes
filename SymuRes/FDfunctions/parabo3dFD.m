function P = parabo3dFD(n,param)
% P = parabo3dFD(n,param)
% Compute the 3D parabolic Fundamental Diagram function
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
%---- P : matrix [P_c; P_b], car and bus productions [veh.m/s]

if ( size(n,1) ) == 1
    n = [n(1,:); zeros(size(n(1,:)))];    
end

vfc = param(1)/param(7);
vfb = param(2)/param(8);
bcc = param(3)/param(7)^2;
bbc = param(4)/param(7)^2;
bcb = param(5)/param(8)^2;
bbb = param(6)/param(8)^2;

P_c = n(1,:).*( vfc + bcc.*n(1,:) + bbc.*n(2,:) );
P_b = n(2,:).*( vfb + bcb.*n(1,:) + bbb.*n(2,:) );

P_c( P_c < 0 ) = 0;
P_b( P_b < 0 ) = 0;

P = [P_c; P_b];

end