function pfit = fitfct3(x,y,f,Npts,pdim,cvgcrit,maxiter)
% pfit = fitfct(x,y,f,pinit,cvgcrit,maxiter)
% Fit the data (x,y) with a function f(x,p) with a set of parameters p
%
% INPUTS
%---- x,y     : vectors of same length, data point set
%---- f       : function, e.g. for a linear fit: f = @(x,p) p(1)*x + p(2)
%               x is the variable and p is the parameter set
%---- pdim    : 2-by-2 matrix, [kmin kmax ; qmin qmax]. Initial value of
%               each parameter is automatically determined thanks to these
%               given ranges for density and flow
%---- cvgcrit : scalar << 1, convergence criterion
%---- maxiter : integer, maximum number of iterations
%
% OUTPUTS
%---- pfit : vector, set of the best parameters found to fit the data

Nparam = 2*Npts + 1;
pinit0 = zeros(1,Nparam);
pdim0 = zeros(Nparam,2);

% pinit0(1) = pinit(1); % kj
% pinit0(2:(Npts+1)) = pinit(2)*ones(1,Npts); % all kc
% pinit0((Npts+2):Nparam) = pinit(3)*ones(1,Npts); % all qc

% pdim0(1,:) = pdim(1,:); % kc
% pdim0(2:(Npts+1),:) = ones(Npts,1)*pdim(2,:); % all kc
% pdim0((Npts+2):Nparam,:) = ones(Npts,1)*pdim(3,:); % all qc

% Chosen range for densities
kmin = pdim(1,1); % minimum acceptable value for all kc and kj
kmax = pdim(1,2); % maximum acceptable value for all kc and kj
% Chosen range for flows
qmin = pdim(2,1); % minimum acceptable value for all qc
qmax = pdim(2,2); % maximum acceptable value for all qc

% The initial values of density are equally distributed over the chosen range
distribk = linspace(kmin,kmax,Npts+1);
pinit0(1) = distribk(Npts+1); % kj, must be greater than all kc
pinit0(2:(Npts+1)) = distribk(1:Npts); % all kc in increasing order

% The initial values of flow are randomly distributed over the chosen range
pinit0((Npts+2):Nparam) = qmin + (qmax - qmin)*rand(1,Npts);

% Acceptable ranges
pdim0(1,:) = [kmin kmax]; % kc
pdim0(2:(Npts+1),:) = ones(Npts,1)*[kmin kmax]; % all kc
pdim0((Npts+2):Nparam,:) = ones(Npts,1)*[qmin qmax]; % all qc

pfit = fitfct2(x,y,f,pinit0,pdim0,cvgcrit,maxiter);

end