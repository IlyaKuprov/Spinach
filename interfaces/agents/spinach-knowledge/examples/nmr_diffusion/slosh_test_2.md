# examples/nmr_diffusion/slosh_test_2.m

- Signature: `slosh_test_2()`

## Purpose

Probability density sloshing around a harmonic oscillator in the presence of a gravitational pull twards the left. Calculation time: seconds.

## Physical / mathematical content

- Diffusion examples. The dominant mathematics is diffusion or advection-diffusion PDE propagation, sometimes with additional spin phase accumulation under gradients.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Probability density sloshing around a harmonic oscillator
- in the presence of a gravitational pull twards the left.
- Calculation time: seconds.
- Set oscillator parameters
- Get the Hamiltonian
- Get the initial state
- Get the propagator
- Run the evolution
