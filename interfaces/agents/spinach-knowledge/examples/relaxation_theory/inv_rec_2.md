# examples/relaxation_theory/inv_rec_2.m

- Signature: `inv_rec_2()`

## Purpose

Simulates inversion-recovery proton spectra for the strychnine spin system at six recovery delays. The source estimates the calculation time in minutes.

## Physical / mathematical content

The spin system is loaded with `strychnine({'1H'})` and set to 14.1 T. Redfield relaxation uses `dibari` equilibrium, kite retention, temperature 298, and a 200 ps correlation time. The script sets a proximity cutoff of 4.0. The basis uses `sphten-liouv`, IK-2 approximation, scalar-coupling connectivity, and proximity level 3; Krylov propagation is disabled.

## Numerical / algorithmic content

For each delay in `[0.01, 0.1, 0.5, 1, 5, 10]` seconds, the equilibrium state is inverted with a 180-degree `Ly` pulse, allowed to recover under the rotating-frame NMR Hamiltonian plus `1i*relaxation`, and tipped by a 90-degree `Ly` pulse. The script acquires a proton FID with sweep width 6500, 8192 points, offset 2800 Hz, and zero-fills to 65536 points. It applies exponential apodisation with parameter 5, Fourier transforms the signal, and plots the real spectrum.

## Implementation structure

The code builds the strychnine spin system and selected basis, sets acquisition parameters, computes `rho_eq` and an `L+` proton detection coil, then forms the rotating-frame Liouvillian. A six-iteration loop applies the inversion, recovery evolution, read pulse, and acquisition; each result is apodised, transformed, and drawn in one panel of a 2-by-3 figure labelled by recovery delay.
