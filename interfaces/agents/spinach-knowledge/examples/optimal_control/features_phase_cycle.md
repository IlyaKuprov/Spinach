# examples/optimal_control/features_phase_cycle.m

- Signature: `features_phase_cycle()`

## Purpose

This example configures a state-transfer pulse from proton `Lz` to transverse fluorine magnetisation and a two-row phase cycle. Its stated design check is that changing the phase of the fluorine channel should produce the corresponding phase change in the final ¹⁹F magnetisation. This describes the intended test; source reading does not establish its outcome.

## Spin model and transfer

The three-spin model contains `1H`, `13C` and `19F` at 9.4 T; each chemical shift is 0.0 ppm. The H–C and C–F scalar couplings are 140 Hz and −160 Hz. The basis is `sphten-liouv` with no approximation. The normalised initial state is `Lz` on proton spin 1; the target is the normalised sum of `L+` and `L-` on fluorine spin 3.

## Pulse design and phase-cycle test

The control set is `Lx/Ly` on all three isotopes, with channel map `[1;1;2;2;3;3]`. The five configured power levels span 0.8–1.2 × 10³ × 2π rad/s. Fifty slices of 0.2 ms give a 10 ms pulse. The script uses the `SNS` penalty with weight 100, `lbfgs`, and a maximum of 100 iterations; the initial guess is a random `6 x 50` array. The phase-cycle matrix has two rows: all-zero phases, then `[0 0 0 pi pi]`. The test code reads column 4 and rotates the fluorine control rows 5:6 by that phase before each simulation.

After optimising, the script simulates each phase-cycle row and calls `stateinfo` on the resulting state. It does not print a scalar fidelity or include observed output, so the intended phase response and convergence are not reported here.

Source: [examples/optimal_control/features_phase_cycle.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/features_phase_cycle.m).
