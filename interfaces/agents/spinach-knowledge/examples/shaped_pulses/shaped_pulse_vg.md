# examples/shaped_pulses/shaped_pulse_vg.m

- Signature: `shaped_pulse_vg()`

## Purpose

Apply a Veshtort-Griffin E1000B 90-degree selective pulse to a system of 31 proton spins with nearest-neighbor J couplings and linear coupling topology. Calculation time: seconds.

## Physical / mathematical content

- The magnetic field is 14.1; the 31 `1H` spins have scalar Zeeman shifts linearly spaced from -4 to 4, with nearest-neighbor scalar couplings of 10.
- The pulse is applied with a 480 Hz frequency offset through `H+2*pi*480*Lz`, starting from an `Lz` state.

## Numerical / algorithmic content

- The basis uses `sphten-liouv` formalism, `IK-2` approximation, `scalar_couplings` connectivity, and proximity level 1.
- `vg_pulse('E1000B',500,duration)` generates amplitudes for a 0.01 s pulse on a 500-step time grid; `shaped_pulse_xy` propagates the state using `expv-pwc`.
- Acquisition uses a 5000 Hz sweep and 2048 points. The FID receives exponential apodisation with parameter 6, followed by an FFT with 16384-point zero filling and `fftshift`.

## Implementation structure

- Create the spin system and basis, apply the `nmr` assumptions, and construct the Hamiltonian and `1H` control and offset operators.
- Execute the shaped pulse, acquire the FID with `liquid` and `@acquire`, then plot the imaginary part of the spectrum in Hz.
