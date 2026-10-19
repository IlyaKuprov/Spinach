# multiclean.m

`multiclean(A,B,tol)` computes a matrix product in output-column panels and applies the existing `clean_up` rounding to each complete panel before retaining it. This limits temporary unpruned fill-in; it does not round partial inner-product sums. It is useful where an ordinary `A*B` followed by cleanup would briefly retain far more entries than the cleaned result.

## Syntax

```matlab
C=multiclean(A,B,tol);
```

`A` and `B` are conformable two-dimensional single or double arrays, sparse or full, on the CPU or GPU, subject to MATLAB’s native multiplication restrictions. `tol` is a finite non-negative real CPU scalar. A zero tolerance uses native multiplication without cleanup. Positive tolerances use the same rounding operation as `clean_up`, not a simple magnitude cutoff.

The helper preserves the native product’s device and sparse/full representation; precision follows `clean_up`, so a single-precision tolerance can lower a double-precision result. It has no spin-system argument and therefore cannot apply custom density thresholds or a cleanup-disable flag. Spinach callers select zero tolerance when cleanup is disabled and retain their final `clean_up(spin_system,C,tol,true)` to apply the whole-matrix storage policy.

Panel sizing depends on input storage and, on the GPU, available device memory. It bounds explicit raw panels, not undocumented multiplication workspace, final result storage, or all assembly temporaries. Runtime can increase. Blocked dense reductions can differ in their last bits from a single native multiplication; values at rounding boundaries need not be bitwise identical.

Do not introduce cleanup between factors, inside sums, or before normalization merely to use this helper: moving a rounding boundary changes the numerical algorithm. Matrix-free objects and polyadics retain their own multiplication paths.

## Related functions

- `clean_up`: rounding and spin-system-dependent storage policy.
- `propagator`, `expmint`: scaling-and-squaring and auxiliary-matrix integrals with potentially large temporary fill-in.
