# examples/esr_sol_pulsed/soft_3_pulse_deer_2e.m

- Signature: `soft_3_pulse_deer_2e()`

## Purpose

Sets up a three-pulse, soft-pulse DEER simulation for two electron spins. The source identifies the soft-pulse treatment as Fokker–Planck formalism and estimates calculation time in minutes.

## Spin system

The two electrons use `sys.magnet=0.3451805`. The first Zeeman tensor is `[2.284, 2.123, 2.075]` with Euler angles `[135, 90, 45]`°; the second is `[2.035, 2.013, 1.975]` with angles `[30, 60, 120]`°. Their coordinates are `[0,0,0]` and `[20,0,0]` Å. The example enables `sodd` (the source labels this as spin-orbit corrections to the DD couplings), uses the `sphten-liouv` formalism without a basis approximation, and disables trajectory-level SSR algorithms.

## Sequence and calculation

The initial state and detection coil are electron `Lz` and `L+`. Powder averaging uses `rep_2ang_6400pts_sph`, propagation method `expm`, and `verbose=0`. The EPR parameters are offset `-4e8`, sweep `3e9`, 256 points, zero-fill 2048, axis units `GHz-labframe`, derivative 0, inverted axis enabled, and assumptions `'deer'`. The source does not annotate units for the numeric offset and sweep values; `GHz-labframe` is the explicit axis-units setting.

The three pulses have ranks `[2,2,2]`, durations `[20,50,40]` ns, phases `[π/2,π/2,π/2]`, power entries `2π × [8e6,8e6,8e6]`, and frequency entries `[9.720e9,10.255e9,9.720e9]`. Pulse power and frequency units are not separately annotated in this file. The first-to-third-pulse gap is 1 μs; the second-pulse scan has 100 steps; the echo time is 100 ns with 100 points.

Execution and diagnostic plotting are delegated to `deer_3p_soft_diag(spin_system,parameters)`. The script itself does not store a returned trace or define figure axes; the result is therefore the diagnostic output of that routine for the configured DEER sequence, not a numerical result embedded in this example.

Source code: [`examples/esr_sol_pulsed/soft_3_pulse_deer_2e.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/soft_3_pulse_deer_2e.m).
