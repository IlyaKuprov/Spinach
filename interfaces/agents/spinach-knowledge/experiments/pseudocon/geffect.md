# experiments/pseudocon/geffect.m

- Signature: `g=geffect(spin_system,states)`

## Purpose

Computes the effective g-tensor for a user-selected Kramers doublet, following [Equations 61 and 62](https://doi.org/10.1063/1.4793736).

## Parameters / inputs

- `spin_system` — Spinach spin system.
- `states` — indices of the selected states, numbered in ascending energy order from the lowest-energy state.

## Outputs

- `g` — 3-by-3 g-tensor matrix in Bohr magneton units.

## Method

The routine obtains each spin's g-tensor and constructs its spin operators and the total magnetic-moment components. It builds and symmetrises the lab-frame Hamiltonian at orientation `[0 0 0]`, diagonalises it, sorts the eigenstates by ascending energy, and selects `states`. The magnetic-moment operators are projected into that selected subspace; the matrix `G` is assembled using Equation 61, and the returned tensor is `real(sqrtm(G))).

## References

- [10.1063/1.4793736](https://doi.org/10.1063/1.4793736)
- [Spin Dynamics Wiki: geffect.m](https://spindynamics.org/wiki/index.php?title=geffect.m)
