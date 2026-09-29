# examples/relaxation_theory/sle_esr_nitroxide_1.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/sle_esr_nitroxide_1.m) · Signature: `sle_esr_nitroxide_1()`

## Purpose and spin model

Compares two calculated ESR spectra for a nitroxide model: one from stochastic Liouville equation (SLE) propagation and one from Bloch–Redfield–Wangsness (BRW/Redfield) relaxation. These are model calculations, not experimental spectra. The source loads `../standard_systems/nitroxide.log` through `gparse`, then calls `g2spinach` with the spin-selection argument `{{'E','E'},{'N','14N'}}`, `[0 0]`, and `options.no_xyz=1`. It sets `sys.magnet=3.5` and a proximity cutoff of `4.0`; the source does not state units for these values. Both calculations use the full `sphten-liouv` basis with no approximation.

## SLE and Redfield calculations

The SLE setup uses `max_rank=10`, `tau_c=5e-11` (unit not stated), `L+` on `E` as both initial state and detection coil, no decoupling, and `E` as the selected spin. Its sweep is `[-3e8,-1e8]`, with `1024` points, `1024` zero-fill points, `GHz-labframe` axis units, and axis inversion enabled. It calls `gridfree(...,'esr')` with `@slowpass`, processes the result with `fdvec(...,5,1)`, and plots the real spectrum.

For the comparison, the source switches to `inter.relaxation={'redfield'}`, zero equilibrium, secular retention, and `inter.tau_c={5e-11}`, then rebuilds the spin system and basis. The BRW path uses `liquid(...,'esr')` with `@slowpass` and the same displayed acquisition settings; it too applies `fdvec(...,5,1)` and plots the real spectrum. The two panels are labelled `SLE` and `BRW`.
