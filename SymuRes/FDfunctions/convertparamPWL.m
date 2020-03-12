function paramPWL = convertparamPWL(param)
% paramPWL = convertparamPWL(param)
% Convert abscissa-ordinate (kc,qc) FD points into slope-intercept (wc,qc0)
% parameters to define a piecewise linear FD
%
% INPUTS
%---- param : Npts-by-2 matrix, points to define the piecewise linear FD
%             column 1 contains the densities kc in increasing order [veh/m]
%             column 2 contains the flows qc corresponding to each density [veh/s]
%
% OUTPUTS
%---- paramPWL : (Npts-1)-by-2 matrix, points to define the piecewise linear FD
%                column 1 contains the slopes wc [m/s]
%                column 2 contains the q-intercepts qc0 corresponding to each slope [veh/s]


kc = param(:,1); % abscissa of branch start point [veh/m]
qc = param(:,2); % ordinate of branch start point [veh/s]

Npts = length(kc);
wc = zeros(1,Npts-1);
qc0 = zeros(1,Npts-1);
for i = 1:(Npts-1)
    wc(i) = (qc(i+1) - qc(i))/(kc(i+1) - kc(i));
    qc0(i) = (qc(i)*kc(i+1) - qc(i+1)*kc(i))/(kc(i+1) - kc(i));
end

paramPWL = [wc; qc0]';

end