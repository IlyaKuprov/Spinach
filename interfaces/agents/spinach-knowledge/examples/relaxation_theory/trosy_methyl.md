# examples/relaxation_theory/trosy_methyl.m

- Signature: `trosy_methyl()`

## Purpose

Simulates methyl TROSY in a rapidly rotating `13CH3` group of a slowly tumbling protein using the Fokker-Planck formalism. The source notes a calculation time of minutes.

## Physical / mathematical content

- The model contains three equally populated methyl rotamers, each with one `13C` and three `1H` spins. Their coordinates and shielding-derived chemical shift tensors are assigned with cyclic permutations of the proton positions.
- Absolute shielding tensors are converted to traceless chemical shift tensors using `-remtrace(shielding{n})`. The proton tensors receive guessed isotropic offsets of `0.8`, `1.0`, and `1.2`.
- Scalar couplings within each rotamer are `125` for each carbon–proton pair and `-12` for each proton–proton pair.
- Methyl turning is represented by a three-state exchange-rate matrix with diagonal entries `-2*k_jump` and off-diagonal entries `k_jump`, where `tau_m=1e-11` and `k_jump=1/(2*tau_m)`.

## Numerical / algorithmic content

- The magnetic field is `14.1`. The calculation uses `sphten-liouv` formalism with no basis approximation and disables `zte` for high accuracy.
- Both frequency-domain calculations use `gridfree(spin_system,@slowpass,parameters,'nmr')`, `tau_c=50e-9`, `max_rank=3`, no decoupling, and matching `L+` initial and detection states for the observed isotope.
- The `13C` spectrum uses a sweep of `[-300 300]` Hz and `1024` points; the `1H` spectrum uses `[200 1000]` Hz and `2048` points. The real spectra are plotted side by side.

## Implementation structure

- Defines Cartesian coordinates and DFT absolute shielding tensors for one carbon and three protons, then constructs the corresponding chemical shift tensors and scalar-coupling matrix.
- Builds a 12-spin system partitioned into three four-spin rotamers, assigns their coordinates, shift tensors, and intrarotamer couplings, and sets equal populations and exchange rates.
- Creates the spin system, selects its basis, calculates the `13C` and `1H` spectra, and plots them with `plot_1d`.
