# examples/fundamentals/nutation_dist_test.m

- MATLAB implementation: [examples/fundamentals/nutation_dist_test.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/nutation_dist_test.m)

- Signature: `nutation_dist_test()`

## Purpose

A synthetic forward-and-inverse demonstration of recovering an RF-field (nutation-frequency) distribution from a nutation curve with same-coil excitation and detection. The source uses reciprocity weighting: each ensemble member's detected transverse signal is multiplied by its own RF amplitude before the ensemble signals are combined. The script generates the curve itself; it does not read an experimental measurement.

## Spin model and synthetic input

The Spinach system contains one `1H` isotope, `sys.magnet = 14.1`, zero scalar Zeeman shift, the `sphten-liouv` formalism, and basis approximation `none`. It sets `Lz` as the initial state, `Lx` and `Ly` as detection channels, and `Lx` as the RF operator; the drift Hamiltonian comes from `hamiltonian(assume(spin_system,'nmr'))`.

The source samples 201 RF angular frequencies from 25 to 65 kHz, represented in rad/s. The normalised bimodal density combines a 0.8-weight Gaussian centred at 50 kHz with 3.0 kHz width and a 0.2-weight Gaussian centred at 38 kHz with 2.5 kHz width. For each grid point, it propagates the spin system under the drift plus RF term, with `dt = 2e-6 s` and `npts = 256`, then accumulates the complex transverse response with probability-mass and RF-amplitude weighting.

The resulting curve is phase-shifted by 1.9 radians and normalised by its maximum magnitude. The source then adds complex Gaussian noise with scale `2e-3` after `rng(1)`, making the noise sequence reproducible for the same MATLAB random-number implementation.

## Distribution recovery and output

The inverse call is `nutation_dist(curve,dt,lambda)` with `lambda = 3e2`; the source identifies this as second-derivative Tikhonov regularisation. It plots the source and recovered probability densities against nutation frequency in kHz. The baseline source comment estimates calculation time in seconds; that comment is not a measured runtime from this review.

## Scope and limitations

The example is a one-spin, on-resonance synthetic recovery with a selected distribution, noise level, and regularisation parameter. It contains no numerical recovery-error metric or pass/fail assertion, and it does not establish performance on experimental data or other distributions. No MATLAB execution or recovered numerical values are claimed here.
