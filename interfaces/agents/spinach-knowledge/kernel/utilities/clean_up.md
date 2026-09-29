# kernel/utilities/clean_up.m

## Purpose

Array clean-up utility. Drops non-zero elements with magnitude below the user-specified tolerance and converts between sparse and full storage depending on the density of non-zeroes in the array.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/clean_up.m>

## Behaviour

- Syntax: `A=clean_up(spin_system,A,nonzero_tol)`.
- Objects of class `opium` are returned unchanged.
- If `nonzero_tol` is `0` or `NaN`, the function returns immediately without cleaning.
- Cell arrays are processed recursively, element by element.
- `polyadic` objects are processed recursively through their `prefix`, `suffix`, and `cores` fields.
- Consistency is enforced by a local `grumble` subfunction: `A` must be numeric, `nonzero_tol` must be numeric, and `nonzero_tol` must be a positive real scalar; violations raise errors.
- Cleaning is skipped entirely when the string `'clean-up'` is present in `spin_system.sys.disable`.
- When enabled, the generic method applies `A=nonzero_tol*round((1/nonzero_tol)*A)`, which snaps values to the nearest multiple of the tolerance and thereby zeroes entries smaller than the tolerance.
- Storage conversion rules:
  - A small sparse matrix with at least one non-zero and any dimension smaller than `spin_system.tols.small_matrix` is converted to full.
  - A big sparse matrix whose non-zero fraction `nnz(A)/numel(A)` exceeds `spin_system.tols.dense_matrix` is converted to full.
  - A big full matrix whose non-zero fraction is below `spin_system.tols.dense_matrix` and whose dimensions all exceed `spin_system.tols.small_matrix` is converted to sparse.

## Inputs and outputs

Inputs:

- `spin_system` — spin system object supplying `sys.disable` and `tols.small_matrix` / `tols.dense_matrix` thresholds.
- `A` — a numerical array or a cell array thereof (also supports `polyadic` objects).
- `nonzero_tol` — nonzero tolerance; a positive real numeric scalar.

Output:

- `A` — cleaned-up array.

## References

- Spinach Wiki page: <https://spindynamics.org/wiki/index.php?title=clean_up.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/clean_up.m>
