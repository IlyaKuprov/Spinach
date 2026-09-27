# examples/shaped_pulses/shaped_pulse_q5.m

- Signature: `shaped_pulse_q5()`

## Purpose

90-degree Q5 pulse on a chain of 31 strongly coupled protons. Calculation time: seconds.

## Physical / mathematical content

- At a magnetic field of 14.1, the system contains 31 `1H` spins with scalar Zeeman values spanning -4 to 4 and nearest-neighbor scalar couplings of 10.
- The pulse acts on an initial `Lz` state using `Lx` and `Ly` controls, with a 480 Hz frequency offset applied through `H+2*pi*480*Lz`.

## Numerical / algorithmic content

- The `sphten-liouv` basis uses the `IK-2` approximation, `scalar_couplings` connectivity, and proximity level 1.
- A 200-point waveform is read from `q5_1000.pk` for a 0.012-second pulse. Its amplitude is calibrated using `8*(pi/2)*npoints/(sum(A)*duration)`, converted from amplitude and phase to Cartesian controls, and propagated with `shaped_pulse_xy` using `expv-pwc`.
- Acquisition uses a 5000 Hz sweep, 2048 points, and zero filling to 16384 points. The FID receives exponential apodisation with parameter 6; the imaginary part of its shifted Fourier spectrum is plotted.

## Implementation structure

- Create the spin system, construct its basis, apply the `nmr` assumptions, and obtain the Hamiltonian and `Lx`, `Ly`, and `Lz` operators.
- Prepare the `Lz` state, execute the calibrated Q5 pulse, then acquire with an `L+` coil, no decoupling, and zero acquisition offset.
- Apodise and Fourier-transform the FID, then plot the resulting spectrum on a Hz axis.
