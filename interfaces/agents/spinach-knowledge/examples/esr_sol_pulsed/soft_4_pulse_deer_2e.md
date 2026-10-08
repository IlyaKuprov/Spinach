# examples/esr_sol_pulsed/soft_4_pulse_deer_2e.m

- Function: `soft_4_pulse_deer_2e()`.

## Model

This example sets up a two-electron, four-pulse DEER simulation. Its header identifies the soft pulses as simulated with the Fokker-Planck formalism and estimates a runtime of minutes. The two isotopes are `{'E','E'}` at `sys.magnet=0.3451805`. The Zeeman principal-value triplets are `[2.284 2.123 2.075]` and `[2.035 2.013 1.975]`; the Euler-angle triplets are `[135 90 45]` and `[30 60 120]` degrees, converted to radians in the assignments. Coordinates are `[0 0 0]` and `[20 0 0]` Angstrom, so the entered separation is 20 Angstrom along x. `sys.enable={'sodd'}` enables the source-commented spin-orbit corrections to dipole-dipole couplings.

## Sequence and acquisition

The basis is `sphten-liouv` with `approximation='none'`; trajectory-level SSR is disabled. Initial state and detection operator are `Lz` and `L+` for `E`. The orientation grid is `rep_2ang_6400pts_sph` and the method is `expm`. EPR controls are `offset=-4e8`, `sweep=3e9`, 256 points, zero-fill 2048, and `axis_units='GHz-labframe'`; derivative is off, `invert_axis=1`, and `assumptions='deer'`. Raw offset and sweep units are not stated in this file.

All pulse ranks are 2; phases are `[pi/2 pi/2 pi/2 pi/2]`. Durations are `[20,40,50,40] ns`, coded as `[20e-9 40e-9 50e-9 40e-9]` seconds. The power vector is `2*pi*[8e6 8e6 8e6 8e6]` and the frequency vector is `[9.720e9 9.720e9 10.255e9 9.720e9]`. The source does not explicitly label units for raw frequency or power values. Echo controls are `p1_p2_gap=0.5e-6` s (0.5 microseconds), `p2_p4_gap=1.5e-6` s (1.5 microseconds), `p3_nsteps=100`, `echo_time=120e-9` s (120 ns), and `echo_npts=240`.

The final call is `deer_4p_soft_diag(spin_system,parameters)` under a “Simulation and plotting” comment. This caller defines the spin system, initial/detection operators, pulses, and echo sampling, but does not itself expose a returned signal or specify plotted axes; those details belong to the called engine. The source contains no numerical output.

## Source

[examples/esr_sol_pulsed/soft_4_pulse_deer_2e.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/soft_4_pulse_deer_2e.m)
