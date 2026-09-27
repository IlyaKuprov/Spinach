# examples/esr_sol_pulsed/soft_3_pulse_deer_3e.m

- Signature: `soft_3_pulse_deer_3e()`

## Purpose

Simulates three-pulse DEER for a three-electron system, with soft pulses treated using the Fokker–Planck formalism. The source estimates a calculation time of minutes.

## Spin system

- Sets the magnetic field to 0.3451805 and uses three electron spins (`'E'`).
- Specifies anisotropic Zeeman principal values and Euler angles for all three electrons.
- Places the electrons at `[0 0 0]`, `[20 0 0]`, and `[0 0 20]` Å, and enables spin-orbit corrections to the dipolar couplings (`'sodd'`).
- Uses the `sphten-liouv` basis without approximation and disables `trajlevel`.

## Sequence and computation

- Uses `Lz` as the initial state and `L+` as the detection state, with the `rep_2ang_6400pts_sph` orientation grid and `expm` propagation.
- Sets three rank-2 pulses with durations `[20 50 40]` ns, phases of π/2, powers of `2π × 8 MHz`, and frequencies `[9.720 10.255 9.720]` GHz.
- Sets a 1 µs first-to-third-pulse gap, 100 steps for the second pulse, and a 50 ns echo time sampled at 100 points.
- Configures a 256-point EPR sweep with 2048-point zero filling, then runs the simulation and plotting through `deer_3p_soft_diag(spin_system,parameters)`.