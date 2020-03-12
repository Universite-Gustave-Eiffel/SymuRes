%% DEFAULT PARAMETERS
%--------------------------------------------------------------------------
% Set default parameters for the Simulation and Assignment structures
% (when not already set in SimulSettings)


% Default parameters - Simulation
%--------------------------------------------------------------------------
if ~isfield(Simulation,'MergeModel')
    Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
end
if ~isfield(Simulation,'DivergeModel')
    Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
end
if ~isfield(Simulation,'TripbasedSimuFactor')
    Simulation.TripbasedSimuFactor = 0.5; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
end
if ~isfield(Simulation,'OutfileDebug')
    Simulation.OutfileDebug = 0; % 0: save essential information / 1: save all simulation information for debug (heavier files)
end


% Default parameters - Assignment
%--------------------------------------------------------------------------
if ~isfield(Assignment,'MinimumGap')
    Assignment.MinimumGap = 0.01; % Threshold for Wardrop equilibrium convergence
end
if ~isfield(Assignment,'NumViolationsThreshold')
    Assignment.NumViolationsThreshold = 0.05; % Threshold of the relative difference between the assignment coefficients
% of two consecutive iterations to consider that a route is in violation
end
if ~isfield(Assignment,'NumViolationsTolerance')
    Assignment.NumViolationsTolerance = 0.05; % Threshold of the proportion of routes that are considered in violation
% (convergence criterion)
end
if ~isfield(Assignment,'MaxIteration')
    Assignment.MaxIteration = 10; % Maximum number of iterations for the convergence loop
end
if ~isfield(Assignment,'WeightMSWA')
    % For more information check the MSWA-I (Liu et al. 2007, DOI:10.1007/s11067-007-9023-x)
    Assignment.WeightMSWA = 2; % weight (> 0) for the MSWA
end
if ~isfield(Assignment,'GammaMSWA')
    % For more information check the MSWA-I (Liu et al. 2007, DOI:10.1007/s11067-007-9023-x)
    Assignment.GammaMSWA = 0; % initialization parameter
end
if ~isfield(Assignment,'ConvergenceIndicatorID')
    % Convergence Indicator:
    % 1 - Gap criterion + Maximum iterations.
    % 2 - Number of violations + Maximum iterations.
    % 3 - Gap criterion + Number of violations + Maximum iterations.
    Assignment.ConvergenceIndicatorID = 3;
end
if ~isfield(Assignment,'model')
    % Assignment models ID:
    % 1   - Deterministic User Equilibrium (DUE)
    % 2   - Stochastic User Equilibrium (SUE)
    % 100 - Manual (user-defined coefficients in Assignment.ManualCoefficients)
    % 101 - Equi-probability
    % 102 - Pro-rata of the number of micro trips
    % 103 - Weighting on reservoirs (weights in Assignment.ReservoirWeight)
    % 104 - Penalty on reservoirs (penalties in Assignment.PenalizedResPath and .PenalizedAssignCoeff)
    Assignment.model = 1;
end
if ~isfield(Assignment,'Behavior')
    % Behavior models ID:
    % 1 - Utility maximizer (rational)
    % 2 - Satisficing behavior (boundedly rational)
    % 3 - Regret minimizer (cf Regret theory)
    Assignment.Behavior = 1;
end
if ~isfield(Assignment,'UtilityID')
    % Utility models ID (for SUE):
    % 1 - Uncertainty on the trip lengths
    % 2 - Uncertainty on the mean speed
    % 3 - Uncertainty on the trip lengths and mean speed
    % 4 - Uncertainty on the trip lengths, mean speed and inhomogeneities on the MFD
    Assignment.UtilityID = 2;
end
if ~isfield(Assignment,'parameters')
    % Parameters definition:
    %
    % DUE: There are no parameters needed to be defined.
    % BR-DUE: Assignment.ALparam has to be defined and corresponds to the
    % Indifference Band. It should be defined between 0 and +Inf.
    %
    %
    % SUE or BR-SUE:
    % Uncertainty on the trip lengths:
    % Assignment.parameters(1) - ID of the utility: Assignment.parameters(1)=1.
    % Assignment.parameters(2) - Number of trip length draws inside
    % each reservoir.
    %
    % Uncertainty on the mean speed:
    % Assignment.parameters(1) - ID of the utility: Assignment.parameters(1)=2.
    % Assignment.parameters(2) - Number of mean speed draws inside
    % each reservoir.
    %
    % Uncertainty on the trip lengths and mean speed:
    % Assignment.parameters(1) - ID of the utility: Assignment.parameters(1)=3.
    % Assignment.parameters(2) - Number of trip length draws inside
    % each reservoir.
    % Assignment.parameters(3) - Number of mean speed draws inside
    % each reservoir.
    % Assignment.parameters(4) - Total number of draws for the Utility function
    % total combination of mxn draws).
    %
    % Uncertainty on the trip lengths, mean speed and inhomogeneities on the MFD:
    % Assignment.parameters(1) - ID of the utility: Assignment.parameters(1)=4.
    %
    Assignment.parameters(1) = 10000;
    Assignment.parameters(2) = 10000;
    Assignment.parameters(3) = 100000;
end
if ~isfield(Assignment,'MinimumGap')
    Assignment.MinimumGap = 0.01; % Threshold for Wardrop equilibrium convergence
end