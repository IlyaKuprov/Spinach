# examples/shaped_pulses/shaped_pulse_chirp_xy.m

- Signature: `shaped_pulse_chirp_xy()`

## Purpose

Demonstrate a chirped inversion pulse and plot its waveform alongside the resulting band-selective inversion spectrum. Calculation time: seconds.

## Physical / mathematical content

- Models 31 `1H` spins at a 14.1 T magnetic field, with scalar Zeeman shifts spanning `-4` to `4` and scalar couplings of `20` between adjacent spins.
- Generates an adaptive WURST chirp with `chirp_pulse(500,0.1,2000,20,'wurst-adaptive')`, returning two transverse waveform components and their segment durations.
- Applies the shaped pulse to the initial `Lz` state, homospoils the result with the `'destroy'` option, then applies a global `pi/2` pulse about `Ly` before acquisition through an `L+` coil.

## Numerical / algorithmic content

- Uses the `sphten-liouv` formalism with `IK-2` approximation, scalar-coupling connectivity, and proximity level `1`.
- Propagates the two-component waveform using `shaped_pulse_xy` with the `'expv-pwc'` method.
- Acquires `2048` points over a `5100` Hz sweep at zero offset, applies exponential apodisation with parameter `6`, then computes the real part of an `8192`-point, `fftshift`-centred FFT.

## Implementation structure

- Constructs the spin system, basis, initial state, detection coil, Hamiltonian, relaxation and kinetics operators, and transverse `Lx`/`Ly` operators.
- Generates and applies the chirp, homospoils, applies the global hard pulse, and acquires and processes the signal.
- Plots both chirp components against cumulative segment duration and displays the processed band-selective inversion spectrum.
