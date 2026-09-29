# examples/nmr_diffusion/slosh_test_2.m

- MATLAB implementation: [examples/nmr_diffusion/slosh_test_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_diffusion/slosh_test_2.m)

- Signature: `slosh_test_2()`

## Purpose

This companion to slosh_test_1 animates one-dimensional wavefunction propagation in a harmonic oscillator with a nonzero gravitational setting. The source describes the pull as leftward and estimates a runtime of seconds. It is a particle-wave example, not an NMR simulation.

## Oscillator and initial state

The oscillator parameters are a force constant of 2.00×10³ N/m, mass 1.00 kg, box size 2.00 m, 100 grid points, and gravitational acceleration 250 m/s². The Hamiltonian and coordinate grid come from `oscillator`; the initial state is `exp(-50*(xgrid-0.1).^2)`, centred at 0.1 on that grid. The script does not explicitly normalise this initial state.

## Propagation and plot

The code forms `expm(full(-1i*H*0.001))` and applies that same propagator 1000 times. Each animation frame plots `abs(psi)+xgrid.^2` in red and `xgrid.^2` in blue, with the horizontal range [−1, 1] m and vertical range [0, 1.5]. The labels identify particle coordinate in metres and probability density in arbitrary units; the red trace is the source's amplitude magnitude plus the reference curve, rather than an explicitly squared-modulus density.

Only an interactive plot is produced. No saved numerical result, normalisation check, or quantitative comparison is included in the script.
