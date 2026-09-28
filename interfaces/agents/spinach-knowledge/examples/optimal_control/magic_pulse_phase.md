# examples/optimal_control/magic_pulse_phase.m

- Signature: `magic_pulse_phase()`

## Purpose

A template for optimising a broadband ¹³C 90-degree excitation pulse that tolerates resonance offsets and RF power calibration errors. The example targets a 28.18 T magnet, approximately 200 ppm (60 kHz) of excitation bandwidth, and negligible effects from a worst-case ¹³C–¹H J-coupling of about 200 Hz. The source comment states a duration limit of 1/(100J) = 50 µs, whereas the implemented pulse has 60 intervals of 1 µs each (60 µs total); the code and comment therefore do not agree. The desired transfers are {Lz → Lx, Ly → Ly, Lx → −Lz}, with RF nutation frequencies from 50 to 70 kHz across the coil. Estimated calculation time: minutes.

- Reference: http://dx.doi.org/10.1016/j.jmr.2005.12.010
- Source contacts: ilya.kuprov@weizmann.ac.il; david.goodwin@inano.au.dk

## Physical / mathematical content

- The ensemble comprises 100 non-interacting ¹³C spins at equally spaced chemical shifts from −100 to +100 ppm. The `sphten-liouv` formalism with `IK-2`, proximity level 1, and `scalar_couplings` connectivity retains complete single-spin bases while ignoring multi-spin orders in this case.
- Normalised `Lx`, `Ly`, and `Lz` states define three simultaneous transfers: `Lx → −Lz`, `Ly → Ly`, and `Lz → Lx`. The control operators are `Lx` and `Ly`; the drift Hamiltonian is obtained under the `nmr` assumption.
- Phase samples are optimised by `fmaxnewton(spin_system,@grape_phase,guess)` with `control.method='lbfgs'` and a 200-iteration limit. The amplitude profile remains fixed at ones; robustness is sampled at ten RF power levels from 50 to 70 kHz, expressed as angular frequencies by multiplication by `2*pi`.

## Numerical / algorithmic content

- Both control operators map to the ¹³C channel through `control.channels=[1; 1]`. The pulse grid is `1e-6*ones(1,60)`, and the initial phase guess is `(pi/5)*randn(1,60)`. Plotting options are `phi_controls`, `xy_controls`, `robustness`, and `spectrogram`.
- After optimisation, `polar2cartesian` converts the phase profile and an amplitude profile of `mean(control.pwr_levels)*control.amplitudes` into Cartesian controls. `shaped_pulse_xy` simulates their action on an initial `Lz` state using `expv-pwc`.
- The resulting state is acquired on ¹³C with an `L+` detection state, no decoupling, zero offset, a 70,000 Hz sweep, 2,048 points, 16,384-point zero filling, a ppm axis, and axis inversion. The FID receives Gaussian apodisation with parameter 10 before its shifted Fourier transform.
- For comparison, a conventional hard pulse is simulated from `Lz` using zero pulse frequency, phase `pi/2`, power `2*pi*60e3`, duration `4.2e-6`, rank 3, and the `expv` method. Its FID receives the same apodisation and Fourier transform. The real spectra are plotted in separate subplots.

## Implementation structure

1. Set the magnetic field and construct the 100-spin chemical-shift ensemble; create the Spinach system and basis.
2. Prepare and normalise the three spin states, then obtain control operators and the drift Hamiltonian.
3. Configure fixed-amplitude, phase-only control and run LBFGS GRAPE optimisation through `fmaxnewton` and `grape_phase`.
4. Convert the optimised pulse to Cartesian controls, simulate it, acquire and process its spectrum, and plot the result.
5. Acquire and process a conventional hard-pulse spectrum for comparison.