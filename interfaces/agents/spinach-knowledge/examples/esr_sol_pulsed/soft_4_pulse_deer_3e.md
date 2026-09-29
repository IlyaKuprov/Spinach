# examples/esr_sol_pulsed/soft_4_pulse_deer_3e.m

- Function: `soft_4_pulse_deer_3e()`.

## Model

This is the three-electron counterpart of the four-pulse DEER example. The source identifies soft-pulse Fokker-Planck simulation and estimates a runtime of hours. The field setting is `sys.magnet=0.3451805`; isotopes are `{'E','E','E'}`. Zeeman principal-value triplets are `[2.284 2.123 2.075]`, `[2.035 2.013 1.975]`, and `[1.935 1.895 1.895]`. Their Euler-angle triplets are `[135 90 45]`, `[30 60 120]`, and `[60 40 20]` degrees (converted to radians in the code). Coordinates in Angstrom are `[0 0 0]`, `[20 0 0]`, and `[0 0 20]`. Spin-orbit corrections to dipole-dipole couplings are enabled by `sys.enable={'sodd'}`.

## Sequence and acquisition

The basis is `sphten-liouv` with no approximation; trajectory-level SSR is disabled. Initial state and detection operator are `Lz` and `L+` on `E`. The grid is `rep_2ang_6400pts_sph` and the method is `expm`. The source sets `offset=-4e8`, `sweep=3e9`, 256 points, zero-fill 2048, `axis_units='GHz-labframe'`, derivative off, axis inversion on, and the `deer` assumption. Raw offset and sweep units are not stated.

Pulse ranks are all 2; phases are `[pi/2 pi/2 pi/2 pi/2]`; durations are `[20,40,50,40] ns` (coded in seconds); powers are `2*pi*[8e6 8e6 8e6 8e6]`; and frequencies are `[9.720e9 9.720e9 10.255e9 9.720e9]`. The file does not annotate raw frequency or power units. Echo controls are `p1_p2_gap=0.5e-6` s (0.5 microseconds), `p2_p4_gap=1.5e-6` s (1.5 microseconds), `p3_nsteps=100`, `echo_time=120e-9` s (120 ns), and `echo_npts=240`.

The script passes these settings to `deer_4p_soft_diag(spin_system,parameters)` under a “Simulation and plotting” comment. It does not itself define plotted axes or numerical results; the called engine owns those details. The hours figure is the source header's estimate, not a newly measured runtime.

## Source

[examples/esr_sol_pulsed/soft_4_pulse_deer_3e.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/soft_4_pulse_deer_3e.m)
