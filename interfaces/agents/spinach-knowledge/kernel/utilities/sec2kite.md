# kernel/utilities/sec2kite.m

## Purpose

Converts a secular relaxation superoperator into the Redfield "kite" form by dropping all non-longitudinal cross-relaxation processes. This is useful when the relaxation superoperator is very large but TROSY-like effects are negligible.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sec2kite.m>

## Behaviour

- Calls `grumble(spin_system,R)` to enforce consistency before processing.
- Records the nonzero count of the input superoperator (`nnz_before`).
- Uses `lin2lm(spin_system.bas.basis)` to compile the index of all longitudinal product states in the basis; longitudinal states are those where the sum of absolute values of the corresponding rows of `M` is zero.
- Converts `R` to XYZ format via `find(R)`.
- Zeros all rates except self-relaxation and longitudinal cross-relaxation terms, keeping entries where both row and column indices are longitudinal states, or where the row and column indices are equal.
- Recomposes the relaxation superoperator with `sparse(rows,cols,vals,length(R),length(R))` and records the new nonzero count (`nnz_after`).
- Reports to the user that non-longitudinal cross-relaxation processes were dropped and prints the reduction in `nnz(R)` from before to after.

Consistency checks performed by the internal `grumble` function:

- Requires `spin_system.bas.formalism` to be `'sphten-liouv'`; otherwise errors with `this function requires sphten-liouv formalism.`
- Requires `R` to be numeric and square; otherwise errors with `R must be a square matrix.`
- Computes the unit state via `unit_state(spin_system)` and errors with `R appears to be thermalised, cannot proceed.` if `norm(R*unit,2)` exceeds `1e-10`.

## Inputs and outputs

Syntax:

```
R = sec2kite(spin_system, R)
```

Inputs:

- `spin_system` — spin system object whose basis must be in the `sphten-liouv` formalism.
- `R` — relaxation superoperator; must be a square numeric matrix and must not be thermalised.

Outputs:

- `R` — relaxation superoperator in Redfield kite form, with non-longitudinal cross-relaxation processes removed.

## References

- Spin Dynamics Wiki page for `sec2kite.m`: <https://spindynamics.org/wiki/index.php?title=sec2kite.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/sec2kite.m>
