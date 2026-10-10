# examples/shaped_pulses/shaped_pulse_chirp_xy.m

- Signature: `shaped_pulse_chirp_xy()`

## Purpose

This example applies a WURST chirp intended for band-selective inversion, then displays its RF waveform and the processed spectrum from a simulated acquisition. The plotted spectrum is a model output; the script does not report an experimental measurement or a convergence study.

## Spin system and pulse

The model is a 31-proton chain at 14.1 T. The scalar Zeeman shifts run linearly from −4 to +4 ppm, and each adjacent pair is assigned a 20 Hz scalar coupling. The basis is the IK-2 restricted state-space basis with scalar-coupling connectivity at proximity level 1.

The call `chirp_pulse(500,0.1,2000,20,'wurst-adaptive')` requests 500 waveform points over 0.1 s, a 2000 Hz sweep bandwidth, and WURST smoothing exponent 20 with adaptive sampling. It returns the two Cartesian RF components Cx and Cy in rad/s and segment durations in seconds; the example plots the components against cumulative duration. The initial density operator is the 1H Lz state. The waveform is propagated with shaped_pulse_xy using the NMR Hamiltonian, Lx/Ly controls, and the expv-pwc piecewise-constant propagation option.

## Preparation and observable

After the chirp, homospoil with the destroy option removes the remaining transverse coherence in the modeled state, and a global Ly hard pulse of π/2 is applied. Acquisition uses the Hamiltonian, relaxation and kinetics operators constructed by the script, with a 5100 Hz sweep, 2048 points, zero filling to 8192 points, and an axis in Hz. The FID receives exponential apodisation with parameter 6; the real part of its zero-filled Fourier transform is plotted. There is no explicit spatial gradient waveform in this example; homospoil is the stated coherence-destruction step.

The shaped_pulse_xy implementation cites DOI [10.1016/j.jmr.2004.08.017](https://doi.org/10.1016/j.jmr.2004.08.017). The example and implementation can be inspected at [shaped_pulse_chirp_xy.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/shaped_pulses/shaped_pulse_chirp_xy.m) and [shaped_pulse_xy.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/shaped_pulse_xy.m).

This is a single-substance, chemistry-free example: `kinetics` supplies a zero generator. Initial states use `state`, and detection uses unweighted `coil_state`.
