# kernel/summaries/summary_chemistry.m

Source: [kernel/summaries/summary_chemistry.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_chemistry.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_chemistry.m)

- Signature: `summary_chemistry(spin_system)`

## Purpose

When `spin_system.chem.parts` contains more than one subsystem, reports each subsystem's number and the spin indices stored in its parts entry. Within that same multi-subsystem branch, it reports reaction and flux data when those fields are present and contain reportable entries.

## Reported values and units

- If `spin_system.chem.rates` exists, the function removes its diagonal before enumerating entries. It prints a table headed `N(from)`, `N(to)`, and `Rate(Hz)`; each value is formatted with a sign and three digits after the decimal in scientific notation (`%+0.3e`). The two indices are printed from the matrix column and row, respectively.
- If `spin_system.chem.flux_rate` exists, its nonzero entries are enumerated with `find`. A point-to-point flux table is printed only when entries are present; its indices are printed from the matrix row and column, respectively. Values are likewise labelled Hz and formatted as signed scientific notation with three digits after the decimal.

Subsystem and spin indices are counts/indices, not physical units. No subsystem or rate table is produced when there is only one chemical part. Rate-table headings are emitted when the rates field exists, even if removing the diagonal leaves no entries.

## Parameters / inputs

- `spin_system` - Spinach spin system structure.

## Output and side effects

Writes subsystem membership and applicable tables through `report.m` to the console or configured output. It checks that `spin_system` is a structure and otherwise does not modify it.
