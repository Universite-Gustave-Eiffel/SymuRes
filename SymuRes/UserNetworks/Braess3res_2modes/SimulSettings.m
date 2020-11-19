%% SIMULATION SETTINGS
%--------------------------------------------------------------------------

Simulation.Duration = 5000; % Simulation duration [s]
Simulation.TimeStep = 10; % Simulation time step [s]

Simulation.MergeModel = 'demprorata'; % 'demprorata', 'demfifo', 'equiproba' or 'endogenous'
Simulation.DivergeModel = 'maxdem'; % 'maxdem', 'decrdem' or 'queuedyn'
Simulation.TripbasedSimuFactor = 0.7; % factor < 1 to scale down the demand level and increase the trip-based solver computation time
Simulation.OutfileDebug = 0; % 0: save essential information / 1: save all simulation information for debug (heavier files)


%% Reservoir definition
%--------------------------------------------------------------------------
% Launch the ResDef module
ResDef

Simulation.MFDfct = MFDfct;
Simulation.Entryfct = Entryfct;
Simulation.Exitfct = Exitfct;


%% Assignment parameters
%--------------------------------------------------------------------------
% Assignment is a structure that gathers all the parameters of the assignment process

% Assignment by fixed periods
Assignment.Periods = [0 Simulation.Duration]; % List of times when traffic assignment is updated [s]
% -- Set to [0 Simulation.Duration] for only one assignment period

% Route calculation
Assignment.PredefRoute = 0; % 1: the routes are already provided (no assignment) / 0: calculate the routes according to the OD matrix
% -- See demand definition in DemDef.m

% Convergence
Assignment.Convergence = 0; % 1: MSA convergence loop for each assignment period / 0: no convergence loop
Assignment.MinimumGap = 0.01; % Threshold for Wardrop equilibrium convergence
Assignment.NumViolationsThreshold = 0.05; % Threshold of the relative difference between the assignment coefficients
% of two consecutive iterations to consider that a route is in violation
Assignment.NumViolationsTolerance = 0.05; % Threshold of the proportion of routes that are considered in violation
% (convergence criterion)
Assignment.MaxIteration = 10; % Maximum number of iterations for the convergence loop
% Parameters for the MSWA. For more information check the MSWA-I (Liu et al. 2007, DOI:10.1007/s11067-007-9023-x)
Assignment.WeightMSWA = 2; % weight (> 0) for the MSWA
Assignment.GammaMSWA = 0; % initialization parameter

% Convergence Indicator:
% 1 - Gap criterion + Maximum iterations.
% 2 - Number of violations + Maximum iterations.
% 3 - Gap criterion + Number of violations + Maximum iterations.
Assignment.ConvergenceIndicatorID = 3;

% Definition of the macroscopic choice set:
Assignment.NumShortestPath = 3; % Max number of shortest Path that will be considered during the assignment process
Assignment.ShortestPath_Threshold = 0.001; % This threshold allows to eliminate the less sampled micro path to define the macro path
Assignment.N_path_samples = 10000;

% Trip lengths model ID:
% Model 0: Assume a Gaussian distribution of trip lengths in reservoirs
% Model 1: Aggregation by reservoir
% Model 2: Aggregation by reservoir and destination
% Model 3: Aggregation by origin, reservoir and destination
% Model 4: Aggregation by macro-path
Assignment.TripLengthMethod = 0; % ID method calculation: 1, 2, 3 or 4
Assignment.LengthStd = 300; % standard deviation of the trip length distribution for method 0 [m]
% In case the user wants to update the trip lengths for method 1 or 2.
Assignment.TripLengthUpdate = 0; % 1 - Yes; 0 - No;

% Assignment models ID:
% 1   - Deterministic User Equilibrium (DUE)
% 2   - Stochastic User Equilibrium (SUE)
% 100 - Manual (user-defined coefficients in Assignment.ManualCoefficients)
% 101 - Equi-probability
% 102 - Pro-rata of the number of micro trips
% 103 - Weighting on reservoirs (weights in Assignment.ReservoirWeight)
% 104 - Penalty on reservoirs (penalties in Assignment.PenalizedResPath and .PenalizedAssignCoeff)
Assignment.model = 1;

% Behavior models ID:
% 1 - Utility maximizer (rational)
% 2 - Satisficing behavior (boundedly rational)
% 3 - Regret minimizer (cf Regret theory)
Assignment.Behavior = 1;

% Utility models ID (for SUE):
% 1 - Uncertainty on the trip lengths
% 2 - Uncertainty on the mean speed
% 3 - Uncertainty on the trip lengths and mean speed
% 4 - Uncertainty on the trip lengths, mean speed and inhomogeneities on the MFD
Assignment.UtilityID = 2;

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

Assignment.ALparam = 0.5; % Relative value for the indifference band.
% Assignment.parameters(1) = Assignment.UtilityID;
Assignment.parameters(1) = 10000;
Assignment.parameters(2) = 10000;
Assignment.parameters(3) = 100000;



%% Demand definition
%--------------------------------------------------------------------------
% Launch the DemDef module

DemDef


