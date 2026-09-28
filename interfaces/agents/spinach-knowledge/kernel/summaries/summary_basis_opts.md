# kernel/summaries/summary_basis_opts.m

- Signature: `summary_basis_opts(spin_system)`

## Purpose

Prints a summary of the basis-set options for a Spinach system. Syntax: `summary_basis_opts(spin_system)`.

## Physical / mathematical content

Reports the selected basis formalism. For `sphten-liouv`, it also reports the configured approximation and its correlation-level, connectivity, proximity-cutoff, and interaction-cutoff settings as applicable.

## Numerical / algorithmic content

The approximation report distinguishes `IK-0`, `IK-1`, `IK-2`, `IK-DNP`, `IK-SBS`, and `none`; the source defines a separate report for each option.

## Parameters / inputs

- `spin_system` - Spinach spin system description object.

## Outputs

- Prints through `report.m` to the console or the user-specified output.

## Implementation structure

- Checks that the input is a structure and reports one of the supported formalisms: Zeeman wavefunction, Zeeman Hilbert-space matrix, Zeeman Liouville-space matrix, or spherical-tensor Liouville-space matrix. For `sphten-liouv`, reports the configured approximation: `IK-0` sets an intra-substance correlation order; `IK-1` and `IK-2` report interaction-graph and proximity criteria; `IK-DNP` reports electron/electron-nuclear/nuclear correlation levels; `IK-SBS` reports boson-boson, spin-boson, and spin-spin levels; `none` reports a complete basis on all spins. Relevant interaction cutoffs and proximity distances are printed from the system settings.
