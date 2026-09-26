# examples/nmr_liquids/noe_strychnine.m

- Signature: `noe_strychnine()`

## Purpose

Inversion-recovery NOE effect spectrum on strychnine, with the rightmost proton signal inverted and a pulse-acquire experiment performed after a 500 ms mixing time. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- Inversion-recovery NOE effect spectrum on strychnine, with the rightmost
- proton signal inverted and a pulse-acquire experiment performed after a
- 500 ms mixing time.
- Calculation time: minutes
- Read the spin system properties
- Magnet field
- Disable Krylov propagation
- Basis set
- Relaxation theory parameters
- Proximity cut-off
- Spinach housekeeping
- Build the relaxation superoperator
