# examples/nmr_diffusion/slosh_test_1.m

- MATLAB implementation: [examples/nmr_diffusion/slosh_test_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_diffusion/slosh_test_1.m)

- Signature: `slosh_test_1()`

## Purpose

This short animation shows a one-dimensional wavefunction evolving in a harmonic-oscillator potential without gravity. It is a particle-wave propagation example, not an NMR experiment; the source estimates a runtime of seconds.

## Oscillator and initial state

The call to `oscillator` builds the Hamiltonian and coordinate grid from a force constant of 2.00×10³ N/m, particle mass 1.00 kg, box size 2.00 m, 100 grid points, and gravitational acceleration 0 m/s². The initial wavefunction is `exp(-50*(xgrid-0.6).^2)`, centred at 0.6 on that grid. The script does not explicitly normalise it.

## Propagation and plot

A fixed propagator is formed as `expm(full(-1i*H*0.001))` and applied successively for 1000 iterations. At each iteration the red curve is `abs(psi)+xgrid.^2`; the blue curve is the `xgrid.^2` reference. The axes are set to [−1, 1] m horizontally and [0, 1.5] vertically, with labels “particle coordinate, m” and “probability density, a.u.” The plotted red quantity is the source's amplitude magnitude plus the reference curve, not a separately computed squared-modulus probability density.

The output is an on-screen animation only. The script gives no saved trajectory, normalisation check, or quantitative comparison with an analytic oscillator solution.
