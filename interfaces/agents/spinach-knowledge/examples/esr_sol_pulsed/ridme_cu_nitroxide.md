# examples/esr_sol_pulsed/ridme_cu_nitroxide.m

- Signature: `ridme_cu_nitroxide()`

## Purpose

Simulates RIDME for a Cu(II)–nitroxide two-electron system at 1.249 T (Q-band). The main calculation uses brute-force time propagation with Liouville-space powder averaging, retaining g-tensor orientation effects on the dipolar coupling. The example includes extended T1/T2 relaxation; its analytical treatment uses only the isotropic parts of the electron g factors. Calculation time: seconds.

## Spin system and relaxation

The copper and nitroxide g-tensor principal values are [2.056, 2.056, 2.205] and [2.009, 2.006, 2.003], respectively, with zero Euler angles. The spins are separated by 43 Å. T1 values are 35 μs (Cu) and 2 ms (nitroxide); T2 values are 1.5 μs and 1.3 μs. Relaxation is retained in the lab frame with zero equilibrium state.

## Simulation and processing

The code uses the full sphten Liouville-space basis without approximation and disables trajectory-level SSR. It starts from the electron `Lz` state and probes spin 2. RIDME evolution uses a 16 ns timestep, 25 and 188 steps in the two evolution dimensions, a 35 μs mixing time, and the `rep_2ang_800pts_sph` grid. The reported trace sums the real and imaginary components of the PxPxPx, PyPyPx, MxMxPx, and MyMyPx phase-cycle pathways. The script plots the real and imaginary trace components against time in microseconds.
