# kernel/summaries/summary_pbc.m

- Signature: `summary_pbc(spin_system,header)`

## Purpose

Prints the periodic-boundary-condition vectors stored in a Spinach spin system.

## Physical / mathematical content

Each vector is reported as its three Cartesian components, taken directly from `spin_system.inter.pbc`.

## Numerical / algorithmic content

The routine prints one row for each cell entry in `spin_system.inter.pbc`; components are formatted to three decimal places. It checks that `spin_system` is a structure and `header` is a character string.

## Parameters / inputs

- spin_system - Spinach spin system description object
- header - a string of text to precede the summary

## Outputs

- Prints the vector table through `report.m` to the console or user-specified output.

## Implementation structure

After validation and table headings, the function traverses the periodic-boundary-condition cell array and reports each vector's X, Y, and Z components.
