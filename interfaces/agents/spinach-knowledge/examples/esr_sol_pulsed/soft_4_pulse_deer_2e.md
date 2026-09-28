# examples/esr_sol_pulsed/soft_4_pulse_deer_2e.m

- Signature: `soft_4_pulse_deer_2e()`

## Purpose

Four-pulse DEER simulation for a two-electron system. Soft pulses are simulated using the Fokker-Planck formalism. Calculation time: minutes.

## Physical / mathematical content

- The system contains two electron spins at a magnetic field of 0.3451805 T, positioned 20 Å apart along the x-axis. Their anisotropic Zeeman interactions have principal values `[2.284 2.123 2.075]` and `[2.035 2.013 1.975]`, with Euler angles `[135 90 45]` and `[30 60 120]` degrees, respectively.
- Spin-orbit corrections to the dipolar couplings are enabled with `sys.enable={'sodd'}`. The initial state is electron `Lz`, and detection uses electron `L+`.

## Numerical / algorithmic content

- The calculation uses the `sphten-liouv` basis without approximation, the `rep_2ang_6400pts_sph` orientation grid, and the `expm` propagation method; trajectory-level output is disabled.
- The EPR settings specify a −400 MHz offset, 3 GHz sweep, 256 points, and zero filling to 2048 points on a `GHz-labframe` axis. The calculation uses `deer` assumptions, no derivative, and an inverted axis.
- Four rank-2 pulses have durations of 20, 40, 50, and 40 ns, each with phase π/2 and power 2π × 8 MHz. Their frequencies are 9.720, 9.720, 10.255, and 9.720 GHz. The first-to-second pulse gap is 0.5 μs, the second-to-fourth pulse gap is 1.5 μs, and the third pulse is varied over 100 steps. Echo acquisition uses a 120 ns echo time and 240 points.

## Implementation structure

- Defines the magnetic field, isotopes, Zeeman tensors, electron coordinates, basis, and algorithmic options; then constructs the spin system with `create` and `basis`.
- Sets sequence, EPR, pulse, and echo-timing parameters, then calls `deer_4p_soft_diag(spin_system,parameters)` for simulation and plotting.
