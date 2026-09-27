# examples/esr_sol_pulsed/oop_eseem_photosystem_1.m

- Signature: `oop_eseem_photosystem_1()`

## Purpose

Powder-averaged, time-domain Liouville-space simulation of two-pulse out-of-phase ESEEM for the spin-correlated electron pair [P700⁺,A1⁻] in Photosystem I. It is intended to reproduce Figure 3a of [doi:10.1021/bi048445d](http://dx.doi.org/10.1021/bi048445d). Calculation time: seconds.

## Physical model

The field is 0.3249 T. The P700⁺ g-tensor eigenvalues (principal axes, angles omitted) are [2.00304, 2.00262, 2.00232], from [doi:10.1016/S0005-2728(01)00198-0](http://dx.doi.org/10.1016/S0005-2728(01)00198-0). The A1⁻ g-tensor eigenvalues are [2.00670, 2.00560, 2.00240], from [doi:10.1007/BF00019589](http://dx.doi.org/10.1007/BF00019589). The two electron spins are separated by the coordinates given in the script, with an isotropic dipolar coupling set from 0.6 mT. A diagonal damping model uses a rate of 1.7×10⁶ s⁻¹ and zero equilibrium state.

## Simulation

The calculation uses the full sphten Liouville-space basis without approximation and disables trajectory-level SSR. It starts from `Lz⊗Lz`, detects the first electron's `L+`, and specifies `Ly` as the pulse operator. The powder-averaged out-of-phase ESEEM signal is sampled at 200 points with a 20 ns timestep on the `rep_2ang_1600pts_sph` grid. The plotted trace is −Im(FID) against time in microseconds.
