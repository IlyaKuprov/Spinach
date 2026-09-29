# examples/optimal_control/mas_powder_control.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/mas_powder_control.m)

## Purpose

This example designs a phase-modulated 87Rb pulse to transfer Lz to Ly in a quadrupolar spin system under magic-angle spinning (MAS). It builds a powder-orientation drift ensemble and evaluates the target-state overlap across that ensemble and configured RF offsets and power levels.

## Spin model and MAS conditions

The model uses 87Rb at 9.413 T and a quadrupolar interaction passed to `eeqq2nqi` as `(1.68e6, 0.2, 3/2, [0 0 0])`. The example uses the `sphten-liouv` formalism without basis approximation. Its MAS settings specify the axis `[sqrt(2/3) 0 sqrt(1/3)]`, a -20 kHz spinning rate (the comment identifies the Bruker direction), the `rep_2ang_100pts_sph` powder grid, maximum spinner rank 8, and rotating-frame order 3 for 87Rb.

## Control design and evaluation

Drift Liouvillians are generated with `singlerot` for the qNMR experiment. Normalised Lz and Ly states are lifted over the classical subspace, defining the transfer objective for each drift. The phase-only GRAPE design keeps the amplitude profile fixed and uses 100 intervals of 0.5 microseconds, which the source describes as one rotor period. Three per-channel RF power levels are configured by `2*pi*[110 120 130]*1e3/sqrt(2)`; five equally spaced offset values span -1,000 to +1,000 (the source does not state their unit). The phase initial guess is constant at pi/2, the optimiser is configured as `lbfgs` through `fmaxnewton` and `@grape_phase`, and the iteration limit is 100.

After pulse design, the script propagates the waveform for each drift Liouvillian and computes the real target-state overlap, then displays the mean fidelity. This is the example's evaluation procedure; no numerical fidelity or convergence result is reported here. The source comments estimate a calculation time of hours. Source contacts: ilya.kuprov@weizmann.ac.il and m.carravetta@soton.ac.uk.
