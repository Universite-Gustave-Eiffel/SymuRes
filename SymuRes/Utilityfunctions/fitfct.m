function pfit = fitfct(x,y,f,pinit,pdim,cvgcrit,maxiter)
% pfit = fitfct(x,y,f,pinit,cvgcrit,maxiter)
% Fit the data (x,y) with a function f(x,p) with a set of parameters p
%
% INPUTS
%---- x,y     : vectors of same length, data point set
%---- f       : function, e.g. for a linear fit: f = @(x,p) p(1)*x + p(2)
%               x is the variable and p is the parameter set
%---- pinit   : vector, set of initial parameters p (all must be non null)
%---- pdim    : matrix, rows are the parameters and columns the min and max values
%---- cvgcrit : scalar << 1, convergence criterion
%---- maxiter : integer, maximum number of iterations
%
% OUTPUTS
%---- pfit : vector, set of the best parameters found to fit the data

Ndata = length(x); % number of data points
Nparam = length(pinit); % number of parameters in f

rmse = sqrt(sum((y - f(x,pinit)).^2)/Ndata); % root-mean-square error (rmse) of the fit
% rmse = sum(abs(y - f(x,pinit)))/Ndata;

pfit = pinit; % set the parameters to the initial set
oldrmse = 0;
niter = 0; % counter for iteration number
while abs(rmse - oldrmse)/oldrmse > cvgcrit && niter < maxiter
    for i = 1:Nparam % loop on all parameters
        newp = pfit;
        newp(i) = min((1 + 0.01)*pfit(i),pdim(i,2)); % increase of 1% of the parameter value
        newrmse = sqrt(sum((y - f(x,newp)).^2)/Ndata); % new rmse with this new parameter value
        %newrmse = sum(abs(y - f(x,newp)))/Ndata;
        if newrmse < rmse % if this rmse is better than the current rmse
            pfit(i) = newp(i); % keep the new value of the parameter in memory
            oldrmse = rmse;
            rmse = newrmse; % update the rmse value
        end
        
        newp = pfit;
        newp(i) = max((1 - 0.01)*pfit(i),pdim(i,1)); % decrease of 1% of the parameter value
        newrmse = sqrt(sum((y - f(x,newp)).^2)/Ndata); % new rmse with this new parameter value
        %newrmse = sum(abs(y - f(x,newp)))/Ndata;
        if newrmse < rmse % if this rmse is better than the current rmse
            pfit(i) = newp(i); % keep the new value of the parameter in memory
            oldrmse = rmse;
            rmse = newrmse; % update the rmse value
        end
    end
    niter = niter + 1;
end

end