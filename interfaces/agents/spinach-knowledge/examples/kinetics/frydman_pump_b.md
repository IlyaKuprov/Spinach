# examples/kinetics/frydman_pump_b.m

- Signature: `frydman_pump_b()`

## Purpose

Lucio Frydman's water exchange based spin-lock pump, Figure 9 from https://doi.org/10.1016/j.jmr.2021.107083 Calculation time: seconds

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Lucio Frydman's water exchange based spin-lock pump, Figure 9
- from https://doi.org/10.1016/j.jmr.2021.107083
- Calculation time: seconds
- Number of water protons
- Magnet field
- Core spin system
- Add water protons
- Chemical shifts
- Scalar couplings
- Estimated relaxation times, seconds
- Relaxation theory
- Basis set
