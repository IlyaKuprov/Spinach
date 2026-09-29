# kernel/summaries/summary_mode_mods.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/summaries/summary_mode_mods.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=summary_mode_mods.m)

## Purpose

Summarise derivatives of spin-spin coupling tensors and effective local fields with respect to dimensionless bosonic-mode displacement coordinates.

## What is reported

Each nonempty derivative entry is listed with its two mode indices, derivative order, target, and norm. Coupling-tensor targets are identified by the affected spin pair; their derivative matrices are measured with the Frobenius norm. Effective-field targets identify the affected spin and use the matrix 2-norm of the derivative vector. Both displayed norms are divided by `2*pi` and reported in Hz. The derivative order is the position `m` in the stored derivative-order list.

## Output and inputs

The routine has no returned output. It sends the supplied header, headings, rows, and final separator through `report` to the console or configured report destination. It checks that `spin_system` is a structure and `header` is a character array; the routine does not assign fields of the input structure.
