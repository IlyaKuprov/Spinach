# experiments/pseudocon/geffect.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/geffect.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=geffect.m) · [Reference DOI 10.1063/1.4793736](https://doi.org/10.1063/1.4793736)

## Purpose

Computes an effective g tensor for a user-selected Kramers doublet, as described by Equations 61 and 62 of [the cited paper](https://doi.org/10.1063/1.4793736). It is a Hilbert-space eigenstate calculation, not a pulse-transfer, DANTE, REDOR, or overtone cross-polarisation simulation.

## Inputs and requirements

- `spin_system` must use the `zeeman-hilb` formalism.
- `states` contains two distinct positive integer state indices. They refer to eigenstates numbered from lowest to highest energy; both indices must be within the compiled Hilbert dimension `bas.offsets(end)`.

For every spin, the routine gets its g tensor with `gtensorof`, builds the `L+`, `L-` and `Lz` operators, and combines them into the three magnetic-moment operators. It obtains the lab-frame Hamiltonian, adds the `[0 0 0]` orientation contribution, symmetrises it, diagonalises and sorts its eigenstates, and projects the magnetic-moment operators into the selected two-state subspace. The 3-by-3 matrix `G` is assembled using Equation 61; the returned tensor is `real(sqrtm(G))`.

## Output

- `g` is a dimensionless 3-by-3 effective g-tensor matrix; the source header’s Bohr-magneton wording does not make this output a magnetic moment.
