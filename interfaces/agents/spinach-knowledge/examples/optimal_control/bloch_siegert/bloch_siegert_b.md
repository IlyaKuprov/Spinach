# examples/optimal_control/bloch_siegert/bloch_siegert_b.m

- Signature: `bloch_siegert_b()`

## Purpose

Bloch-Siegert shift compensation demo for a universal rotation pulse over a range of resonance offsets. It compares optimization with and without BSS correction as control power varies. Calculation time: minutes.

## Physical / mathematical content

- The model uses 100 non-interacting (^{13}mathrm{C}) spins with equally spaced offsets from -100 to +100 ppm at `sys.magnet=28.18`. The `IK-2` basis retains the complete basis on each spin while neglecting multi-spin orders. The desired rotation maps (S_x,S_y,S_z) to (-S_z,S_y,S_x). The comparison illustrates BSS-related fidelity loss when the shift is not included in pulse design.

## Numerical / algorithmic content

- L-BFGS uses 500 iterations maximum and `tol_x=1e-4`; GRAPE-XY pulses have 50 slices. Twenty control powers span (10^{-3}) to 1 times the carbon Zeeman frequency. At each power, pulses are optimized with BSS off and on from a shared random guess, then both are evaluated with BSS enabled. The plotted quantity is terminal infidelity versus relative control power.

## Implementation structure

- Set the field, 100 isotope and offset entries, and the `IK-2` basis with `prox_level=1` and `scalar_couplings` connectivity; construct the system and basis; build normalized spin states, control operators, and drift Hamiltonian; configure the ensemble-independent optimizer; sweep powers, design and evaluate both pulses, and plot the infidelity curves.
