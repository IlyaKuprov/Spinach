# examples/dnp_liq/jdnp/fig_3_time_dep_top_row.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/jdnp/fig_3_time_dep_top_row.m`
- Signature: `fig_3_time_dep_top_row()`
- Total lines: 87

## Purpose

A demonstration that the JDNP effect vanishes when the second electron is removed from the system. Proton polarisation as a function of time for specific external fields is plotted. See also the "bot row" simulation where both electrons are active and the JDNP enhancement is present. Further details in: Calculation time: seconds, line-by-line plotting

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- A demonstration that the JDNP effect vanishes when the second
- electron is removed from the system. Proton polarisation as a
- function of time for specific external fields is plotted. See
- also the "bot row" simulation where both electrons are active
- and the JDNP enhancement is present. Further details in:
- Calculation time: seconds, line-by-line plotting
- Load the spin system
- Kill the second electron
- Experiment parameters
- Magnetic field grid, Tesla
- Get a figure going
- Loop over the fields

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `system_specification()`, `rmfield()`, `kfigure()`, `scale_figure()`, `field_grid()`, `g2freq()`, `create()`, `basis()`, `operator()`, `equilibrium()`, `hamiltonian()`, `assume()`, `relaxation()`, `state()`, `evolution()`, `subplot()`.
