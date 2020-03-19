# SymuRes
A Multi-Reservoir MFD-Based Traffic Simulation Platform

Version 1.0 - 2019

License
-------
SymuRes is licensed under the terms of the GNU GPL v3 license.

Acknowledgments
---------------
This simulation platform was developed during the PhD of Guilhem Mariotte. It also received contributions from the PhD thesis of Sérgio F. A. Batista for the traffic assignment part. It is one of the outcomes of the [MAGnUM project](https://magnum-erc.weebly.com/), held by Prof. Ludovic Leclercq and funded by the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation program (grant agreement No 646592).

Authors
-------
- Guilhem Mariotte - guilhem.mariotte@univ-eiffel.fr | [guilhemmariotte.com](http://guilhemmariotte.com/en)
(Simulation platform design, traffic flow solvers, pre-processing and post-processing modules)

- Sergio F. A. Batista - sab21@nyu.edu
(DTA module, assignment and convergence loop)

Designing a simulation
----------------------
In your network folder in the `UserNetworks/` directory, create or edit the scripts `SimulSettings.m` (global settings of the simulation, i.e. duration, timestep, assignment, etc), `ResDef.m` (characteristics of the reservoirs), `DemDef.m` (definition of the macroscopic OD matrix).

Launching a simulation
----------------------
Specify the network, the solver and the name of the output file in the script `Main.m`. Run this script for a classical simulation (with eventually DTA) in the main SymuRes directory. The output of the simulation is always saved into the `outputs/` folder in the corresponding network folder.

Plotting the results of a simulation
------------------------------------
Specify the simulation outputs to load at the beginning of the script `scrPlotResults.m`. Several simulations can be loaded at the same time. Run parts of this script to plot the desired results (accumulations, flows, N-curves, travel times, etc). For each plot script, indicate the relevant list of simulation results, reservoirs, routes to compare.
