# examples/esr_sol_pulsed/soft_3_pulse_deer_2e.m

- Signature: `soft_3_pulse_deer_2e()`

## Purpose

Three-pulse DEER simulation for a two-electron system using soft pulses and the Fokker–Planck formalism. Calculation time: minutes.

## Spin system

The field is 0.3451805 T. The electron g-tensor principal values are [2.284, 2.123, 2.075] and [2.035, 2.013, 1.975], with Euler angles [135°, 90°, 45°] and [30°, 60°, 120°]. Spin-orbit corrections to the dipolar coupling are enabled. The spins are separated by 20 Å. The calculation uses the full sphten Liouville-space basis without approximation and disables trajectory-level SSR.

## Pulse sequence and acquisition

The initial state and detection operator are the first electron's `Lz` and `L+`. The three pulses have ranks [2, 2, 2], durations [20, 50, 40] ns, phases π/2, and powers 2π×8 MHz. Their carrier frequencies are [9.720, 10.255, 9.720] GHz. The p1–p3 gap is 1 μs; the second-pulse interval uses 100 steps. The echo time is 100 ns with 100 echo points.

The calculation uses the `rep_2ang_6400pts_sph` powder grid and the `expm` method. Acquisition parameters are offset −4×10⁸, sweep 3×10⁹, 256 points, and zero-fill to 2048; the axis is in GHz lab-frame units and inverted. The sequence is run by `deer_3p_soft_diag` with the DEER approximation flag enabled.
