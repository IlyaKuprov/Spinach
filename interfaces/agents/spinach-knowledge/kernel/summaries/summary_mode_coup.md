# kernel/summaries/summary_mode_coup.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_mode_coup.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_mode_coup.m)

## Purpose

Print a table of bosonic-mode couplings stored in a Spinach system.

## What is reported

For each nonempty coupling entry, the table gives the component indices, a coupling-type label, and the stored amplitude divided by `2*pi` and formatted in Hz. The reported categories are exchange, cross-Kerr, longitudinal, radiation pressure, and dispersive. Longitudinal entries are labelled radiation pressure when both endpoint component types belong to C, V, or T; otherwise they retain the longitudinal label. Every nonempty entry in the exchange array is printed as `exchange`; the routine does not give spin-mode exchange a separate label. The source description includes spin-mode exchange among the couplings covered by this summary.

## Output and inputs

The routine has no returned output. It sends the supplied header, table headings, rows, and separators through `report` to the console or configured report destination. The input checks require `spin_system` to be a structure and `header` to be a character array; no input-structure fields are assigned by the routine.
