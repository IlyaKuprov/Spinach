# examples/nmr_diffusion/slosh_test_1.m

- Signature: `slosh_test_1()`

## Purpose

Illustrates time evolution of a wavefunction in a one-dimensional harmonic oscillator with zero gravitational acceleration. The source estimates seconds of calculation time.

## Physical and numerical content

The oscillator is configured with force constant 2e3 N/m, mass 1 kg, a 2 m box, and 100 grid points. The initial state is exp(−50(xgrid−0.6)²). The example obtains the Hamiltonian and coordinate grid from `oscillator`, forms the fixed-step propagator `expm(−1i*H*0.001)`, and applies it 1000 times. At each step it plots abs(psi)+xgrid² alongside the xgrid² reference curve.
