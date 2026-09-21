# examples/dnp_liq/jdnp/fig_1_exch_and_field_scan.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/jdnp/fig_1_exch_and_field_scan.m`
- Signature: `fig_1_exch_and_field_scan()`
- Total lines: 96

## Purpose

Matching condition plot for JDNP -proton polarisation at a particular time as a function of the external field and the inter-electron exchange coupling. Further details in: Calculation time: hours, line-by-line plotting

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Matching condition plot for JDNP -proton polarisation at
- a particular time as a function of the external field and
- the inter-electron exchange coupling. Further details in:
- Calculation time: hours, line-by-line plotting
- Load the spin system
- Experiment parameters
- Field and coupling grids
- Preallocate output array
- Create and scale the figure
- Loop over the fields
- Set the magnet field
- Trityl and free electron frequencies

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `system_specification()`, `kfigure()`, `scale_figure()`, `field_grid()`, `g2freq()`, `exch_grid()`, `create()`, `basis()`, `operator()`, `equilibrium()`, `state()`, `hamiltonian()`, `assume()`, `relaxation()`, `evolution()`, `dnp()`.
