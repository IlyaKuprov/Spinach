# kernel/summaries/summary_rlx_flux.m

`summary_rlx_flux(spin_system)` reports validated point-to-point magnetisation flux rates, their direction, and whether correlations are retained. It works before basis construction and prints nothing for absent or zero flux. Rows of the flux-rate matrix are sources, columns are destinations. `kinetics(spin_system,'report')` calls it before compiling reactions. The function does not construct or change a generator.
