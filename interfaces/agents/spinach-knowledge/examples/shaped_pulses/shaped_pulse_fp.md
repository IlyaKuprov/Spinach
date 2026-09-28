# examples/shaped_pulses/shaped_pulse_fp.m

- Signature: `shaped_pulse_fp()`

## Purpose

Simulates an off-resonance rectangular soft pulse using the Fokker-Planck formalism. The pulse frequency offset accumulates as additional phase during the pulse, as it would during a chirp. Calculation time: seconds.

## Physical / mathematical content

- The system contains 31 `1H` spins at a magnetic field of 14.1 T. Their scalar Zeeman interactions span `-4` to `4`, and adjacent spins have scalar couplings of `10`.
- The background NMR Hamiltonian, `1H` control operators `Lx` and `Ly`, and an initial `Lz` state are passed to `shaped_pulse_af`. The pulse uses a frequency offset of `1922.4` Hz.

## Numerical / algorithmic content

- The basis uses the `sphten-liouv` formalism with `IK-2` approximation, `scalar_couplings` connectivity, and proximity level `1`.
- The pulse call is `shaped_pulse_af(spin_system,H,Lx,Ly,rho,1922.4,50.0,5e-3,-pi/2,2)`.
- Acquisition detects `1H` with an `L+` coil state, zero offset, a `7000` Hz sweep, and `2048` points. The FID is phase-adjusted by `exp(-1i*0.67)`, exponentially apodised with parameter `6`, Fourier transformed with zero filling to `8192` points, and plotted as the real spectrum on a Hz axis.

## Implementation structure

- Creates the spin system and basis, constructs the Hamiltonian and control operators, and applies the soft pulse to the initial state.
- Runs liquid-state acquisition, processes the FID, and plots the resulting spectrum.
