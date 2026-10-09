# examples/shaped_pulses/shaped_pulse_vg.m

- Signature: `shaped_pulse_vg()`

## Purpose

Apply a Veshtort–Griffin E1000B 90-degree selective pulse to a 31-proton chain with nearest-neighbor couplings, and inspect the simulated spectrum under an explicit frequency offset.

## Physical and numerical model

The model has 31 proton spins (1H) at a 14.1 T field, scalar Zeeman values linearly spaced from -4 to 4, and scalar couplings of 10 between adjacent spins. It uses the `sphten-liouv` formalism, `IK-2` approximation, scalar-coupling connectivity and proximity level 1. The initial state is proton longitudinal (Lz) magnetisation.

The pulse duration is 10 ms, divided into 500 equal time intervals. The source obtains a 500-sample amplitude waveform from `vg_pulse('E1000B', 500, duration)`. It applies that waveform on the x control channel with a 480 Hz offset represented by adding `2*pi*480*Lz` to the Hamiltonian; propagation is requested with `expv-pwc`. The source describes this as a 90-degree selective pulse, but does not report a measured or independently validated excitation profile.

## Acquisition and observable

Liquid-state proton acquisition uses L+ as the coil operator, zero acquisition offset, a 5000 Hz sweep, 2048 points and 16384 zero-fill points. The FID is exponentially apodised with parameter 6 and Fourier transformed; the plotted observable is the imaginary part of the spectrum.

## Scope

The frequency offset is part of the pulse Hamiltonian, not a gradient. The example does not specify a relaxation superoperator, gradient encoding, or experimental validation.

The source comment gives a calculation time of seconds.

Source: [examples/shaped_pulses/shaped_pulse_vg.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/shaped_pulses/shaped_pulse_vg.m)
