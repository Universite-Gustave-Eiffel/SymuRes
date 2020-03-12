# SymuRes
A Multi-Reservoir MFD-Based Traffic Simulation Platform

Version 1.0 - 2019

Licence
-------
SymuRes is licenced under the terms of the GNU GPL v3 licence.

Authors
-------
- Guilhem Mariotte - guilhem.mariotte@univ-eiffel.fr
(Simulation platform design, traffic flow solvers, pre-processing and post-processing modules)

- Sergio Batista - 
(DTA module, assignment and convergence loop)

Designing a simulation
----------------------
Edit the scripts SimulSettings.m (global settings of the simulation, i.e. duration, timestep, assignment, etc), ResDef.m (characteristics of the reservoirs), DemDef.m (definition of the macroscopic OD matrix)

Launching a simulation
----------------------
Specify the network, the solver and the name of the output file in the script Main.m. Run this script for a classical simulation (with eventually DTA) in the main folder. The output of the simulation is always saved into the 'outputs/' folder in the corresponding network folder.

Plotting the results of a simulation
------------------------------------
Specify the simulation outputs to load at the beginning of the script scrPlotResults.m. Several simulations can be loaded at the same time. Run parts of this script to plot the desired results (accumulations, flows, N-curves, travel times, etc). For each plot script, indicate the relevant list of simulation results, reservoirs, routes to compare.
