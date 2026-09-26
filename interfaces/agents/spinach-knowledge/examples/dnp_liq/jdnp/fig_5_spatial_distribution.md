# examples/dnp_liq/jdnp/fig_5_spatial_distribution.m

- Signature: `fig_5_spatial_distribution()`

## Purpose

An illustration of the fact that JDNP effect does not vanish on position and orientation averaging in liquid phase. The simula- tion shows proton polarisation at 20 ms for different proton lo- cations around a radical pair with inter-electron exchange cou- pling chosen to acheive the JDNP effect. Details in: Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- An illustration of the fact that JDNP effect does not vanish on
- position and orientation averaging in liquid phase. The simula-
- tion shows proton polarisation at 20 ms for different proton lo-
- cations around a radical pair with inter-electron exchange cou-
- pling chosen to acheive the JDNP effect. Details in:
- Calculation time: seconds
- Load the spin system
- Set magnet field
- Experiment parameters
- Set microwave offset frequency
- Set the exchange coupling
- Specify coordinate arrays
