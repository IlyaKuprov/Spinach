# examples/esr_sol_pulsed/soft_4_pulse_deer_3e.m

- Signature: `soft_4_pulse_deer_3e()`

## Purpose

Simulates four-pulse DEER for a three-electron system with soft pulses using the Fokker–Planck formalism. The source notes a calculation time of hours.

## Model and numerical setup

- Uses three electron spins at coordinates (in Å) `[0,0,0]`, `[20,0,0]`, and `[0,0,20]`, with specified anisotropic Zeeman tensors and Euler angles. The magnetic field is `0.3451805`; spin-orbit corrections to dipole–dipole couplings are enabled with `sys.enable={'sodd'}`.
- Builds a `sphten-liouv` basis with no approximation and disables `trajlevel`. The initial state is `Lz` and the detection state is `L+` for electron spins.
- Uses the `rep_2ang_6400pts_sph` grid and matrix-exponential propagation. EPR settings include a `-4e8` offset, `3e9` sweep, 256 points, and 2048-point zero filling, with a `GHz-labframe` axis.

## Pulse sequence and output

- Specifies four rank-2 pulses with durations `[20,40,50,40]` ns, phases of `π/2`, and angular-frequency powers of `2π × 8 MHz`. Pulse frequencies are `[9.720,9.720,10.255,9.720]` GHz.
- Sets the first-to-second pulse gap to `0.5` µs, the second-to-fourth pulse gap to `1.5` µs, 100 steps for the third pulse, and a 120 ns echo window sampled at 240 points.
- Runs the simulation and plotting through `deer_4p_soft_diag(spin_system,parameters)`.