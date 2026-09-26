# examples/nmr_liquids/noe_four_spin.m

- Signature: `noe_four_spin()`

## Purpose

Inversion-recovery NOE effect spectrum on a simple four-spin system, with the rightmost proton signal inverted and a pulse-acquire experiment per- formed after a very long (five seconds) mixing time. Sequential NOE hops with alternating signs are clearly visible in the result. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Inversion-recovery NOE effect spectrum on a simple four-spin system, with
- the rightmost proton signal inverted and a pulse-acquire experiment per-
- formed after a very long (five seconds) mixing time. Sequential NOE hops
- with alternating signs are clearly visible in the result.
- Calculation time: seconds
- Set the spin system
- Magnet field
- Basis set
- Relaxation theory parameters
- Spinach housekeeping
- Build the relaxation superoperator
- Get thermal equilibrium state
