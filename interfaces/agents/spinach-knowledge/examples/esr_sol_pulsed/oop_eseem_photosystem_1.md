# examples/esr_sol_pulsed/oop_eseem_photosystem_1.m

- Signature: `oop_eseem_photosystem_1()`
- Source: [`examples/esr_sol_pulsed/oop_eseem_photosystem_1.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/oop_eseem_photosystem_1.m)
- Sequence: [`experiments/esr_dipolar/oopeseem.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/oopeseem.m)
- Figure 3a reference: [doi:10.1021/bi048445d](https://doi.org/10.1021/bi048445d)
- P700⁺ g-tensor eigenvalues: [doi:10.1016/S0005-2728(01)00198-0](https://doi.org/10.1016/S0005-2728(01)00198-0); A1⁻ g-tensor eigenvalues: [doi:10.1007/BF00019589](https://doi.org/10.1007/BF00019589)

## Physical model

This is a powder-averaged, time-domain Liouville-space simulation of two-pulse out-of-phase ESEEM for the spin-correlated ([P700^+, A1^-]) pair in Photosystem I, intended to reproduce Figure 3a of the cited paper. It contains two electron spins, with field 0.3249 T. Their principal g values are [2.00304, 2.00262, 2.00232] for P700⁺ and [2.00670, 2.00560, 2.00240] for A1⁻; both tensors have zero Euler angles in this example. The spins are placed at [0,0,0] and [25.35,0,0] in the Spinach coordinate convention (Å). The scalar coupling is set by `mt2hz(0.6e-3)`, converting a 0.6 mT coupling input to Hz. Relaxation is diagonal damping at (1.7×10^6) s⁻¹ with zero equilibrium. The full sphten Liouville-space basis is used without approximation; trajectory-level SSR is disabled.

## ESEEM protocol and observable

The initial state is `Lz⊗Lz`; the detected state is `L+` on the first electron, and `Ly` is the pulse operator. The `oopeseem` helper applies its π/4 pulse and computes the refocused echo trajectory. Powder averaging uses `rep_2ang_1600pts_sph`; 200 points are calculated at a 20 ns timestep. The plotted signal is `−Im(FID)` against the code-defined time `i × timestep/2` (in microseconds), where (i=0,ldots,199); the plot labels intensity in arbitrary units. The script creates a figure but does not save a data file. Its header estimates a calculation time of seconds.
