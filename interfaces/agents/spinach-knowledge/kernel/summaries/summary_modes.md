# kernel/summaries/summary_modes.m

- Signature: `summary_modes(spin_system,header)`

## Purpose

Prints a parameter summary for bosonic modes represented in a Spinach spin system.

## Physical / mathematical content

The table includes modes whose component type is C, V, or T, labelled as cavity, phonon, or transmon, respectively. It reports each mode's isotope, multiplicity, frequency, anharmonicity, damping, and dephasing values.

## Numerical / algorithmic content

The listed frequency, anharmonicity, damping, and dephasing fields are divided by `2*pi` and printed in Hz. The function checks that `spin_system` is a structure and `header` is a character string.

## Parameters / inputs

- spin_system - Spinach spin system description object
- header - a string of text to precede the summary

## Outputs

- Prints the mode table through `report.m` to the console or user-specified output.

## Implementation structure

After validation and table headings, the function iterates over the spin-system components, selects types C, V, and T, maps those type codes to labels, and prints the stored mode parameters.
