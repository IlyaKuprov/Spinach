# kernel/utilities/validate_sym.m

- Signature: `validate_sym(spin_system,bas)`

Validates declared spin-label permutation symmetry against Zeeman tensors, giant-spin coefficients, and couplings. Interactions must match in the laboratory frame; rotations are not applied.

`spin_system` supplies processed interaction data, spin counts, isotopes, and the interaction cutoff. `bas.sym_group` contains group names; matching entries of `bas.sym_spins` contain spin indices. If `sym_group` is present, both fields must be cell arrays of equal length; group names must be character strings and spin-index vectors must contain at least two valid indices.

If symmetry is disabled or no group is declared, validation is skipped. Otherwise, spins in each group must have the same isotope. Each group operation is checked with tolerance `2*pi*spin_system.tols.inter_cutoff` rad/s; a mismatch raises an error. The function returns nothing.

Source: <https://spindynamics.org/wiki/index.php?title=validate_sym.m>