# examples/optimal_control/distortions/distortions_figure_4_top.m

- Signature: `distortions_figure_4_top()`

## Purpose

Figure 4 (top) from the paper by Rasulov and Kuprov:

## Physical / mathematical content

- The example designs a 13C XY pulse for 100 non-interacting spins spanning −100 to +100 ppm at a magnetic field of 28.18. Its initial states are Sx, Sy, and Sz, with targets −Sz, Sy, and Sx.
- The pulse has 125 intervals of 1 μs; the final five are frozen as dead time. Optimisation uses LBFGS-GRAPE with Lx and Ly controls, five RF power levels from 50 to 70 kHz, NS and SNS penalties, and an RLC distortion ensemble with quality factors from 560 to 640.
- The optimised pulse is benchmarked over quality factors from 200 to 1000 and RF nutation frequencies from 40 to 80 kHz. After applying the modeled distortion, the code simulates the pulse and plots the log of infidelity calculated from the three target-state overlaps.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Figure 4 (top) from the paper by Rasulov and Kuprov:
- Set the magnetic field
- Put 100 non-interacting spins at equal intervals
- within the [-100,+100] ppm chemical shift range
- Select a basis set -IK-2 keeps complete basis on each
- spin in this case, but ignores multi-spin orders
- Run Spinach housekeeping
- Set up spin states
- Get the control operators
- Get the drift Hamiltonian
- Define control parameters
- Last 5 slices are dead time
