# examples/shaped_pulses/shaped_pulse_slr.m

- Signature: `shaped_pulse_slr()`

## Purpose

Shinnar-Le Roux band-selective 90-degree excitation pulse on a chain of 31 strongly coupled protons. Calculation time: seconds.

## Physical / mathematical content

- The system contains 31 `1H` spins at a magnetic field of `14.1`, with scalar Zeeman values spanning `-4` to `4` and nearest-neighbor scalar couplings of `10`.
- The initial state is proton `Lz` magnetization. The SLR pulse acts through the proton `Lx` and `Ly` control operators; the acquired signal uses a proton `L+` coil state.

## Numerical / algorithmic content

- Uses the `sphten-liouv` basis with `IK-2` approximation, `scalar_couplings` connectivity, and proximity level `1`.
- Generates a 90-degree excitation waveform with `slr_pulse(256,15e-3,32,pi/2,0.01,0.01)` and applies its `Cx` and `Cy` components using `shaped_pulse_xy` with the `expv-pwc` method.
- Acquires 2048 points with a sweep of 5000 Hz and zero-fills to 16384 points. Applies exponential apodisation with parameter `6`, then computes `fftshift(fft(fid,parameters.zerofill))`.

## Implementation structure

- Creates the spin system, basis, Hamiltonian, control operators, and initial state.
- Generates and applies the excitation SLR pulse, then runs liquid-state NMR acquisition.
- Plots the pulse components against cumulative duration and the magnitude spectrum as band-selective excitation.
