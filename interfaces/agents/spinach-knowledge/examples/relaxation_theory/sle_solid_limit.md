# examples/relaxation_theory/sle_solid_limit.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/sle_solid_limit.m)

## Purpose

This example compares calculated slow-motion ESR line shapes across a set of SLE ranks and correlation times as the stochastic motion slows toward the solid limit. The source lists a calculation time of hours; it does not report experimental data or numerical line positions.

## Spin model and SLE setup

The spin system contains 14N and an electron; the source sets `sys.magnet=0.3343`. The 14N-electron hyperfine tensor is passed through `gauss2mhz` from the matrix `[10.00 3.544 11.70; 3.544 18.00 5.072; 11.70 5.072 30.00]`. The electron Zeeman matrix is `[2.0065794 -0.0007548 -0.0032848; -0.0007548 2.0056940 -0.0006008; -0.0032848 -0.0006008 2.0048920]`. The source gives no other units for these matrix entries. The basis is `sphten-liouv` with no approximation.

The initial state and detection coil are both the electron raising state `L+`; the selected spin is `E`, and the decoupling list is empty. The sweep is `[-2.2e8, 2e8]` on the `GHz-labframe` axis, with 240 points and 240 zero-fill points. The axis is inverted and derivative mode is off. Each trace is generated with `gridfree`, `slowpass`, and ESR mode.

## Rank and correlation-time series

The script pairs four maximum ranks with four correlation times in order: ranks 3, 7, 15, and 30 with correlation-time settings `1e-9`, `1e-8`, `1e-7`, and `1e-6`, respectively. It calculates and plots the real signal for each pair in its own panel. Each panel is titled with its correlation time and labelled with electron Zeeman frequency in GHz. The source defines this comparison; it does not establish convergence or claim specific simulated intensities.
