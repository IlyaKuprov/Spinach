# examples/relaxation_theory/sle_nmr_dd_csa.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/sle_nmr_dd_csa.m)

## Purpose

This example compares a stochastic Liouville Equation (SLE) NMR calculation with a Bloch-Redfield-Wangsness (BRW) calculation for DD-CSA cross-correlation in a 15N-1H protein amide-bond model. Both traces are calculated; the source does not provide a measured spectrum. Calculation time is listed as seconds.

## Spin model and relaxation pathway

The model has 15N and 1H spins; the source sets `sys.magnet=14.1`. The supplied 15N Zeeman matrix is `[14.89 34.35 0; 34.35 145.74 0; 0 0 150.51]`; the 1H matrix is `[30.75 1.31 0; 1.31 22.65 0; 0 0 11.80]`. The scalar-coupling entry is 10, and the two coordinates are `[-0.451455 -0.678015 0]` and `[-1.475290 -0.641823 0]`. The source does not state units for these matrix entries, scalar-coupling value, or coordinates. In this geometry the intended relaxation interference is between the 15N-1H dipole-dipole interaction and 15N chemical-shift anisotropy (DD-CSA), as identified in the source comment. The proximity cutoff is 4.0. The basis is `sphten-liouv` with no approximation.

## SLE and BRW calculations

For SLE, the source sets maximum rank 10 and `parameters.tau_c=5e-9`. The initial state and coil are both `L+` on 15N; 15N is the selected spin and the decoupling list is empty. `gridfree` with `slowpass` runs in NMR mode.

The BRW comparison sets Redfield relaxation, zero equilibrium, secular relaxation retention, and `inter.tau_c={5e-9}`. It uses the same 15N initial state and coil. The SLE spectrum uses a sweep of `[-0.633e4, -0.629e4]` Hz and 2048 points; the BRW calculation uses the same sweep, point count, and Hz axis. The BRW signal is calculated with `liquid` and `slowpass` in NMR mode.

## Plot

The script places the real SLE and BRW signals in side-by-side panels labelled SLE and BRW, with amplitude in arbitrary units shown for the SLE panel. These are the outputs of the two models, not experimental data.
