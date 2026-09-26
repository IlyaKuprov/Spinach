# examples/relaxation_theory/sat_rec_1.m

- Signature: `sat_rec_1()`

## Purpose

A simple saturation-recovery experiment. Calculation time: seconds.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- A simple saturation-recovery experiment.
- Calculation time: seconds.
- Spin system
- Zeeman interactions
- Complete basis set
- Relaxation theory
- Spinach housekeeping
- Initial state -unit
- Detection state
- Static Hamiltonian superoperator
- Pulse operator
- Relaxation superoperator
