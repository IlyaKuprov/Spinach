# examples/dnp_liq/jdnp/fig_2_tau_and_field_traject.m

- Signature: `fig_2_tau_and_field_traject()`

## Purpose

Time evolution plot for JDNP: proton polarisation as a function of time for specific external fields and rotational correlation times. The inter-electron exchange coupling is set to the match- ing condition at each field. Further details in: Calculation time: seconds, line-by-line plotting

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Time evolution plot for JDNP: proton polarisation as a function
- of time for specific external fields and rotational correlation
- times. The inter-electron exchange coupling is set to the match-
- ing condition at each field. Further details in:
- Calculation time: seconds, line-by-line plotting
- Load the spin system
- Experiment parameters
- Magnetic field grid, Tesla
- Correlation time grid, seconds
- Get a figure going
- Loop over the field grid
- Set magnet field
