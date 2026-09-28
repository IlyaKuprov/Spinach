# kernel/utilities/banner.m

- Signature: `banner(spin_system,identifier)`

## Purpose

Internal kernel helper that prints a named console banner through Spinach's reporting mechanism. User calls are discouraged.

## Parameters / inputs

- `spin_system` - Spinach spin-system structure passed to `report`.
- `identifier` - character string selecting one of the supported banners: `version_banner`, `spin_system_banner`, `basis_banner`, `sequence_banner`, `optimcon`, or `optimisation`.

The version banner prints the Spinach v2.13 label and links to the author list, documentation, and book; the other identifiers print section headings for the spin system, basis set, pulse sequence, optimal control, or optimisation. An unknown identifier raises an error.
