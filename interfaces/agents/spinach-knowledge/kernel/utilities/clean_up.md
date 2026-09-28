# kernel/utilities/clean_up.m

- Signature: `A=clean_up(spin_system,A,nonzero_tol)`

## Purpose

Rounds numerical array entries to increments of a tolerance and selects sparse or full storage according to the array's dimensions and density.

## Numerical / algorithmic content

The function returns `opium` objects unchanged and skips cleanup when `nonzero_tol` is zero or NaN. It recursively processes cells and the prefix and suffix elements of polyadic objects. For ordinary arrays, it rounds with `nonzero_tol*round(A/nonzero_tol)`; unless cleanup is disabled in `spin_system.sys.disable`, it then applies the configured `small_matrix` and `dense_matrix` thresholds to choose sparse or full storage.

## Parameters / inputs

- `spin_system` — spin system descriptor containing cleanup settings and storage thresholds
- `A` — numeric array, cell array, or polyadic object; `opium` objects are returned unchanged
- `nonzero_tol` — positive real rounding tolerance; zero or NaN disables cleanup

## Outputs

- `A` — cleaned-up array or object
