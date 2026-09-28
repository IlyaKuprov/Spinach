# examples/optimal_control/steady_orbit/solid_effect_chirp.m

- Signature: `solid_effect_chirp()`

## Purpose

Optimises the phase of a stroboscopic steady-state DNP pulse, using timing and power settings matching the XiX experiment. The source notes that the calculation can take days on a large parallel cluster.

## Model and optimisation

- Models an electron and a proton at 3.35316 T and 80 K, separated by 3.500 Å. It uses the specified Zeeman tensors, distance- and orientation-dependent proton longitudinal relaxation, and a full spherical-tensor Liouville-space basis.
- Uses a powder grid (`rep_2ang_800pts_sph`) and a 94.0 GHz transmitter. The objective is proton longitudinal magnetisation relative to thermal equilibrium, evaluated in steady state.
- Optimises electron-control phase across 720 pulse samples of 0.5 ns each. The following 20 samples form a frozen 0.5 ns-per-sample ringdown, followed by a frozen 167 µs delay. Microwave power levels are `2*pi*linspace(5,25,20)*1e6` rad/s, with offsets of −2, −1, 0, +1, and +2 MHz.
- Uses `rbfgs` with up to 10,000 iterations, a budget of 500, and 240 processes. A 16-tap filter loaded from `hiper_kernel_trans.mat` is normalised to unit absolute DC gain and applied as control distortion.
- Starts from a smoothed 50 MHz chirp over 360 ns, with its phase negated and shifted by 140 MHz. It calls `fmaxnewton` with `grape_phase` and requests robustness and spectrogram plots.