# kernel/summaries/summary_modes.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_modes.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_modes.m)

## Purpose

Print the stored parameters of bosonic components in a Spinach system.

## What is reported

The table includes only components whose type code is C, V, or T, labelled cavity, phonon, or transmon respectively. For each included component it prints the component index, the isotope string from `spin_system.comp.isotopes`, the type label, `spin_system.comp.mults`, and the stored frequency, anharmonicity, damping, and dephasing values. The four parameter values are divided by `2*pi` and printed in Hz. The isotope string is displayed as component metadata; this routine does not calculate a nuclear spin quantum number or isotope-dependent gyromagnetic ratio.

## Output and inputs

The routine has no returned output. It sends the supplied header, headings, rows, and separators through `report` to the console or configured report destination. Its input checks require `spin_system` to be a structure and `header` to be a character array; the routine does not assign fields of the input structure.
