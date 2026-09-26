# examples/dnp_liq/jdnp/fig_3_time_dep_bot_row.m

- Signature: `fig_3_time_dep_bot_row()`

## Purpose

Time evolution plot for JDNP: proton polarisation as a function of time for specific external fields. The inter-electron exchan- ge coupling is set to the matching condition at each field. See also the "top row" simulation where one of the electrons is re- moved to demostrate that JDNP vanishes. Further details in: Calculation time: seconds, line-by-line plotting

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Time evolution plot for JDNP: proton polarisation as a function
- of time for specific external fields. The inter-electron exchan-
- ge coupling is set to the matching condition at each field. See
- also the "top row" simulation where one of the electrons is re-
- moved to demostrate that JDNP vanishes. Further details in:
- Calculation time: seconds, line-by-line plotting
- Load the spin system
- Experiment parameters
- Magnetic field grid, Tesla
- Get a figure going
- Loop over the fields
- Set magnet field
