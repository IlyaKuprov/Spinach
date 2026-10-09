# examples/esr_sol_pulsed/soft_3_pulse_deer_3e.m

- Signature: `soft_3_pulse_deer_3e()`

## Purpose

Sets up a three-pulse, soft-pulse DEER simulation for a three-electron model. The source identifies the soft-pulse treatment as Fokker–Planck formalism and estimates calculation time in minutes.

## Spin system

The three electrons use `sys.magnet=0.3451805`. Their Zeeman principal values are `[2.284, 2.123, 2.075]`, `[2.035, 2.013, 1.975]`, and `[1.935, 1.895, 1.895]`; the corresponding Euler angles are `[135, 90, 45]`°, `[30, 60, 120]`°, and `[60, 40, 20]`°. Coordinates are `[0,0,0]`, `[20,0,0]`, and `[0,0,20]` Å. The example enables `sodd` (identified in the source as spin-orbit corrections to DD couplings), uses the `sphten-liouv` formalism without a basis approximation, and disables trajectory-level SSR algorithms.

## Sequence and calculation

The initial state and detection coil are electron `Lz` and `L+`. Powder averaging uses `rep_2ang_6400pts_sph`, propagation method `expm`, and `verbose=0`. The EPR parameters are offset `-4e8`, sweep `3e9`, 256 points, zero-fill 2048, axis units `GHz-labframe`, derivative 0, inverted axis enabled, and assumptions `'deer'`. The source does not annotate units for the numeric offset and sweep values; `GHz-labframe` is the explicit axis-units setting.

The three pulses have ranks `[2,2,2]`, durations `[20,50,40]` ns, phases `[π/2,π/2,π/2]`, power entries `2π × [8e6,8e6,8e6]`, and frequency entries `[9.720e9,10.255e9,9.720e9]`. Pulse power and frequency units are not separately annotated in this file. The first-to-third-pulse gap is 1 μs; the second-pulse scan has 100 steps; the echo time is 50 ns with 100 points.

Execution and diagnostic plotting are delegated to `deer_3p_soft_diag(spin_system,parameters)`. The script itself does not store a returned trace or define figure axes; the result is therefore the diagnostic output of that routine for the configured DEER sequence, not a numerical result embedded in this example.

## Relation to the two-electron variant

This version adds an electron at `[0,0,20]` Å with its own Zeeman tensor and orientation to the two-electron geometry. The pulse, powder-grid, propagation, EPR, and timing parameters otherwise match the two-electron example; the changed echo time is 50 ns rather than 100 ns.

Source code: [`examples/esr_sol_pulsed/soft_3_pulse_deer_3e.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/soft_3_pulse_deer_3e.m).
