# kernel/summaries/summary_coordinates.m

Source: [kernel/summaries/summary_coordinates.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_coordinates.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_coordinates.m)

- Signature: `summary_coordinates(spin_system,header)`

## Purpose

Prints a coordinate table for all spins, preceded by the caller-provided header. Each row contains the spin index, isotope string from `spin_system.comp.isotopes`, the three stored components from `spin_system.inter.coordinates`, and the spin label from `spin_system.comp.labels`.

## Reported values and units

Coordinate components are formatted with `%+5.3f` and printed under the `X`, `Y`, and `Z` headings. This function does not specify coordinate units, so none are asserted here. Spin indices are dimensionless indices; isotope and label columns are text.

## Parameters / inputs

- `spin_system` - Spinach spin system structure.
- `header` - Character array printed before the table.

## Output and side effects

Writes the header and coordinate table through `report.m` to the console or configured output. It checks that `spin_system` is a structure and `header` is a character array; it does not otherwise modify either input.
