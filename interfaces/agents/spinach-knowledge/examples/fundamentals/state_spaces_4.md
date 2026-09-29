# examples/fundamentals/state_spaces_4.m

[Source code](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/state_spaces_4.m)

## Purpose

This example sets up a powder magic-angle-spinning (MAS) trajectory calculation for isotopically labelled glycine, initially with proton L+ magnetisation and proton L+ detection; `parameters.spins={'13C'}` separately selects the rotor context’s working channel. It is useful for seeing how a full spherical-tensor Liouville-space basis, a finite orientation grid, rotor-phase averaging, and correlation-order analysis are combined in one simulation. The source describes the calculation as taking hours and notes that a GPU can make it faster.

## System and state-space setup

The system is read from ../standard_systems/glycine.log using gparse and g2spinach, with isotope labels 1H, 13C, and 15N and the three numeric arguments [31.5 182.1 264.5] passed to g2spinach. The script sets the magnet field to 14.1 and uses sphten-liouv with approximation='none'. Its longitudinal subspace specification is {{'15N','13C'}}, and its maximum coherence rank is 17. The Krylov tolerance field is set to 1000; the source gives no interpretation of this value.

The source configures a spinning rate of 2000, rotor axis [1 1 1], sweep 1e5, 64 points, and offset 15000. It selects the 13C spin channel and the orientation grid named rep_2ang_200pts_sph. Both the initial state and receiver are proton L+. GPU enablement appears only as a commented-out line, so the script as written does not explicitly enable GPU execution.

## Calculation and analysis

singlerot is called with the traject callback and NMR mode to produce a trajectory. The script then passes `2*max_rank+1` (35 for this setup) to fpl2rho for rotor-phase averaging. Finally, trajan analyzes the trajectory using correlation_order; the displayed plot uses a logarithmic y-axis restricted to 1e-5 through 0.05.

## Output and scope

The function creates a correlation-order figure; it does not return or save a numerical result. The source defines the system, basis, grid name, propagation call, averaging call, and plot limits, but contains no plotted values or reported trajectory features. It therefore supports describing the calculation setup, not asserting a quantitative MAS response or a conclusion about which state-space contributions dominate. The orientation-grid identifier is retained as given; the source does not explain its quadrature construction or accuracy.
