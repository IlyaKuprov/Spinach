# kernel/summaries/summary_basis_opts.m

Source: [kernel/summaries/summary_basis_opts.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_basis_opts.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_basis_opts.m)

- Signature: `summary_basis_opts(spin_system)`

## Purpose

Reports the selected basis formalism. Only for the `sphten-liouv` formalism, it also reports the configured approximation and its associated levels, graph settings, and cutoffs.

## Reported settings

Formalism labels are: `zeeman-wavef` (Zeeman wavefunction), `zeeman-hilb` (Zeeman Hilbert-space matrix), `zeeman-liouv` (Zeeman Liouville-space matrix), and `sphten-liouv` (spherical-tensor Liouville-space).

For `sphten-liouv`, the approximation report is selected by `spin_system.bas.approximation`:

- `IK-0`: reports `inter_level`, described as correlations of all spins up to that order within each chemical substance.
- `IK-1`: reports `inter_level` for correlations between directly coupled spins; `prox_level` for correlations among spins within `prox_cutoff` Angstrom; and the connectivity graph plus `inter_cutoff` in Hz, below which interaction tensors are dropped by the stated criterion (2-norm).
- `IK-2`: reports nearest-neighbour correlations on the coupling graph, the same proximity-level/distance fields, and the connectivity graph with the interaction-tensor 2-norm cutoff in Hz.
- `IK-DNP`: reports the three entries of `inter_level` as maximum inter-electron, electron-nuclear, and inter-nuclear correlation levels; it also reports nearest neighbours on the coupling graph and the interaction-tensor 2-norm cutoff in Hz.
- `IK-SBS`: reports the three entries of `inter_level` as maximum boson-boson, spin-boson, and spin-spin correlation levels, plus the connectivity graph and interaction-tensor 2-norm cutoff in Hz.
- `none`: reports that the starting basis is complete on all spins.

The code prints configured values; it does not calculate or validate the physical meaning of the selected approximation. Proximity is printed in Angstrom and the interaction cutoff in Hz. Correlation levels are dimensionless integers.

## Parameters / inputs

- `spin_system` - Spinach spin system structure.

## Output and side effects

Writes the selected formalism and, when applicable, approximation details through `report.m` to the console or configured output. It checks that `spin_system` is a structure and otherwise does not modify it. Unknown formalism or approximation values raise an error.
