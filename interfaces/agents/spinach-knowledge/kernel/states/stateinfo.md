# kernel/states/stateinfo.m

- Signature: `stateinfo(spin_system,rho,npops)`

## Purpose

Prints a summary of a state vector's norm and its largest-magnitude basis coefficients. It accepts only the `sphten-liouv` formalism, where `spin_system.bas.basis{n}` supplies the spherical-tensor basis labels.

## Inputs

- `rho` is a numeric column vector.
- `npops` is a positive real integer no greater than the vector length.

## Printed values

The function uses `report` to print:

1. The vector's 2-norm, computed as `norm(rho,2)`.
2. The `npops` entries with the largest `abs(rho)`, sorted in descending order of magnitude. Each row names its hosting substance and global spin indices and contains one `(L,M)` label per local spin, the corresponding coefficient, and its one-based position in the basis vector. A spin with `L=0` is printed as `....`.

The coefficient is the vector entry itself, not its squared magnitude; it is formatted with `%+5.3e`. The function reports no physical unit; values retain the numerical scale and any units of the supplied vector. The labels follow the local spin columns of the hosting block, and the final number is the basis-vector index, not a population.

## Side effects and return

This function reports to the console through `report` and has no return value. It does not alter `rho` or the basis.

- Source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/states/stateinfo.m
- Wiki: https://spindynamics.org/wiki/index.php?title=stateinfo.m
- Related: [basis](../basis.md)
