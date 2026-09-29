# kernel/summaries/summary_couplings.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_couplings.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_couplings.m)

## Purpose

Print a table of selected spin-spin coupling tensors from a Spinach system.

## What is reported

A pair is included when the matrix spectral norm satisfies `norm(J,2) > 2*pi*spin_system.tols.inter_cutoff`. The table identifies the two spin indices and prints all three matrix rows, the isotropic value `trace(J)/3`, and the rank-1 and rank-2 component norms. The tensor decomposition uses `mat2sphten` and `sphten2mat`; each displayed matrix row, isotropic value, and rank norm is divided by `2*pi`, so the displayed coupling values are in Hz. The rank norms use the matrix 2-norm.

## Output and inputs

The routine has no returned output. It sends the supplied header, table headings, rows, and separators to `report`, which routes text to the console or the configured report destination. Its input checks require `spin_system` to be a structure and `header` to be a character array; the routine does not update fields of the input structure.
