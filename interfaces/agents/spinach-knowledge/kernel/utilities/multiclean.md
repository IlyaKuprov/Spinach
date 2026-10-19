# multiclean.m

`multiclean(spin_system,A,B,tol)` computes a matrix product in output-column panels and applies the existing `clean_up` rounding to each complete panel before retaining it. This limits temporary unpruned fill-in; it does not round partial inner-product sums. It is useful where an ordinary `A*B` followed by cleanup would briefly retain far more entries than the cleaned result.

## Syntax

```matlab
C=multiclean(spin_system,A,B,tol);
```

`A` and `B` are conformable two-dimensional single or double arrays, sparse or full, on the CPU or GPU, subject to MATLAB’s native multiplication restrictions. `tol` is a non-negative real numeric CPU scalar, or NaN. Integer and infinite tolerances retain native-product `clean_up` arithmetic. Zero or NaN tolerance, or the `clean-up` disable switch in `spin_system.sys.disable`, uses native multiplication without cleanup. Positive tolerances use the same rounding operation as `clean_up`, not a simple magnitude cutoff.

The helper preserves the native product’s device; precision follows `clean_up`, so a single-precision tolerance can lower a double-precision result. `spin_system` supplies the disable switch and global sparse/full storage thresholds. The final whole-matrix representation follows those thresholds, without a second numerical rounding pass. Callers pass their original tolerance directly and do not need a cleanup-bypass block or an outer `clean_up` call.

Panel sizing depends on input storage and, on the GPU, available device memory after the fixed left operand and any dense output have been allocated. It bounds explicit raw panels, not undocumented multiplication workspace, final result storage, or all assembly temporaries. Runtime can increase. Blocked dense reductions can differ in their last bits from a single native multiplication; values at rounding boundaries need not be bitwise identical.

Do not introduce cleanup between factors, inside sums, or before normalization merely to use this helper: moving a rounding boundary changes the numerical algorithm. Custom objects retain their native multiplication and `clean_up` dispatch inside `multiclean`; callers do not need operand-type branches or to materialise those objects.

## Related functions

- `clean_up`: rounding and spin-system-dependent storage policy.
- `propagator`, `expmint`: scaling-and-squaring and auxiliary-matrix integrals with potentially large temporary fill-in.
