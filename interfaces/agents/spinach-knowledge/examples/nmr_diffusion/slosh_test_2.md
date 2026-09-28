# examples/nmr_diffusion/slosh_test_2.m

- Signature: `slosh_test_2()`

## Purpose

Shows a wavefunction evolving in a harmonic oscillator with a nonzero gravitational term (250 m/s² in the source). The source estimates seconds of calculation time.

## Physical and numerical content

The oscillator uses force constant 2e3 N/m, mass 1 kg, a 2 m box, and 100 grid points. Its initial state is exp(−50(xgrid−0.1)²). After obtaining the Hamiltonian and grid from `oscillator`, the example forms `expm(−1i*H*0.001)` and applies this propagator 1000 times. The animation plots abs(psi)+xgrid² and the xgrid² reference curve at each iteration.
