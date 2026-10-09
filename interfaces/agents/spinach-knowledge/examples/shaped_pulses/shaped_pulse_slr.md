# examples/shaped_pulses/shaped_pulse_slr.m

- Signature: `shaped_pulse_slr()`

## Purpose

Build and apply a Shinnar–Le Roux (SLR) 90-degree band-selective excitation pulse in a strongly coupled 31-proton chain, then inspect the simulated proton spectrum.

## Physical and numerical model

The spin system uses 31 proton spins (1H) at a 14.1 T field, scalar Zeeman values spaced from -4 to 4, and scalar couplings of 10 between adjacent spins. The basis is the `sphten-liouv` formalism with the `IK-2` approximation and scalar-coupling connectivity with proximity level 1. The initial density operator is proton longitudinal (Lz) magnetisation.

The waveform generator is called as `slr_pulse(256, 15e-3, 32, pi/2, 0.01, 0.01)`. It returns x- and y-channel controls `Cx` and `Cy`, with durations `durs`; the plotted control amplitudes are labelled in rad/s and cumulative time in seconds. The example applies both quadratures with `shaped_pulse_xy` and the `expv-pwc` method. The source identifies the pulse as a 90-degree excitation; its plotted frequency response is the simulated result, not an experimental validation.

## Acquisition and observable

Liquid-state acquisition uses proton L+ as the coil operator, zero offset, a 5000 Hz sweep, 2048 acquired points and 16384 zero-fill points. The FID receives exponential apodisation with parameter 6 before a zero-filled Fourier transform; the example plots the pulse waveform and the magnitude spectrum labelled as band-selective excitation.

## Scope

This example specifies an excitation waveform and a simulated spectrum. It does not specify a gradient, relaxation model, experimental comparison, or convergence study.

The source comment gives a calculation time of seconds.

Source: [examples/shaped_pulses/shaped_pulse_slr.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/shaped_pulses/shaped_pulse_slr.m)
