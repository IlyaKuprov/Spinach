# examples/nmr_liquids/noesy_sucrose.m

- Signature: `noesy_sucrose()`

## Purpose

NOESY spectrum of sucrose (magnetic parameters computed with DFT). Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- NOESY spectrum of sucrose (magnetic parameters computed with DFT).
- Calculation time: minutes
- Spin system properties (vacuum DFT calculation)
- Magnet field
- Basis set
- Relaxation theory parameters
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- F2 Fourier transform
