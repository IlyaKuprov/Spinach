# examples/relaxation_theory/inv_rec_2.m

- Signature: `inv_rec_2()`

## Purpose

An example of inversion recovery experiment simulation for a strychnine spin system. Calculation time: minutes.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- An example of inversion recovery experiment simulation
- for a strychnine spin system.
- Calculation time: minutes.
- Read the spin system properties
- Magnet field
- Disable Krylov propagation
- Basis set
- Relaxation theory parameters
- Proximity cut-off
- Spinach housekeeping
- Aquisition parameters
- Set up different recovery delays
