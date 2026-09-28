# kernel/summaries/summary_coordinates.m

- Signature: `summary_coordinates(spin_system,header)`

## Purpose

Prints a table of the atomic coordinates stored in the spin-system structure.

## Physical / mathematical content

For each spin, the table reports the three coordinate components from `spin_system.inter.coordinates`, together with the spin index, isotope, and label. The function reports stored coordinates; it does not compute them.

## Numerical / algorithmic content

After checking that `spin_system` is a structure and `header` is a character string, the function sends the heading and table through `report`. It iterates over `spin_system.comp.nspins`; coordinate components are formatted to three decimal places.

## Parameters / inputs

- spin_system - Spinach spin system description object
- header - a string of text to precede the summary

## Outputs

- Prints to the console or to the user-specified output via `report.m`.

## Implementation structure

The routine validates its two inputs, prints the column headings, then emits one row per spin using its isotope, three stored coordinate components, and label.
