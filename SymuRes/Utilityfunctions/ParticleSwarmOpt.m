function pfit = ParticleSwarmOpt(x,y,fct,pinit,pdim,cvgcrit,maxiter,Npart)
% pfit = ParticleSwarmOpt(x,y,fct,pinit,pdim,cvgcrit,maxiter,Npart)
% Fit the data (x,y) with a function f(x,p) with a set of parameters p
% Use the Swarm Particle algorithm
%
% INPUTS
%---- x,y     : vectors of same length, data point set
%---- fct     : function, e.g. for a linear fit: f = @(x,p) p(1)*x + p(2)
%               x is the variable and p is the parameter set
%---- pinit   : vector, set of initial parameters p
%---- pdim    : matrix, rows are the parameters and columns the min and max values
%---- cvgcrit : scalar << 1, convergence criterion
%---- maxiter : integer, maximum number of iterations
%---- Npart   : integer, number of particles which "travel" through the space of parameters
%
% OUTPUTS
%---- pfit : vector, set of the best parameters found to fit the data

Ndata = length(x); % number of data points
Nparam = length(pinit); % number of parameters in f
ymean = mean(y); % mean value of y, for normalization

% Initialization
position = zeros(Npart,Nparam);
velocity = zeros(Npart,Nparam);
chi = 0.729; % required constant in the velocity update calculation (empirical value)

firsttime = 1; % boolean to indicate the first loop
bestparticleID = 1; % particle index with the best parameter position (best fit)

% Initialize the particle positions and velocities (random choice in the parameter space)
for j = 1:Nparam
    for i = 1:Npart
        position(i,j) = pdim(j,1) + rand*(pdim(j,2) - pdim(j,1));
        %position(i,j) = pinit(j)*(0.5 + 1*rand);
        velocity(i,j) = 0.5*position(i,j);
    end
end
bestparticle_rmse = 1000*ones(1,Npart);
bestposition = position;
newfit = 0;

fprintf('\n')

% Start loop
niter = 0;
ncvg = 0; % marker for convergence (stop the loop after 5 consecutive convergence checkout)
while bestparticle_rmse(bestparticleID) > cvgcrit && niter < maxiter && ncvg < 5
    r1 = rand;
    r2 = rand;
    
    for i = 1:Npart % loop on all particles
        rmse = 1/ymean*sqrt(sum((y - fct(x,position(i,:))).^2)/Ndata); % root-mean-square error of the fit, normalized by the mean
        if firsttime == 1
            bestparticle_rmse(i) = rmse;
        end
        
        if rmse < bestparticle_rmse(i)
            bestparticle_rmse(i) = rmse;
            for j =  1:Nparam
                bestposition(i,j) = position(i,j);
            end
            if bestparticle_rmse(i) < bestparticle_rmse(bestparticleID)
                bestparticleID = i;
            end
        end
        
        for j = 1:Nparam
            % velocity adapts (slow down or increase) as the particle
            % position becomes more and more close to the best position
            velocity(i,j) = chi*(velocity(i,j) + 2.1*r1*(bestposition(i,j) - position(i,j)) + 2*r2*(bestposition(bestparticleID,j) - position(i,j)));
            newpos = position(i,j) + velocity(i,j);
            if pdim(j,1) <= newpos && newpos <= pdim(j,2) % if the new position is between the parameter min and max
                position(i,j) = position(i,j) + velocity(i,j);
            end
        end
    end
    
    % convergence checkout
    oldfit = newfit;
    newfit = sum(position);
    if abs(newfit - oldfit)/oldfit < 0.01
        ncvg = ncvg + 1;
    else
        ncvg = 0;
    end
    
    firsttime = 0;
    niter = niter + 1;
    fprintf('%s%d \t %s%f\n','iter=',niter,'RMSE=',bestparticle_rmse(bestparticleID))
end

pfit = bestposition(bestparticleID,:);

end
        
        