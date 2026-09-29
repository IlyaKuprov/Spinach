# kernel/utilities/cheap_norm.m

## Purpose

Returns the cheapest available matrix norm depending on how the matrix is represented: the infinity-norm for GPU arrays, the 1-norm for CPU arrays, and a lower-bound 1-norm estimate for polyadic objects, which can only multiply vectors. The header notes that some norms are vastly more expensive than others and that this function uses the cheapest ones available.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/cheap_norm.m>

## Norm guarantee by representation

For a GPU array, `cheap_norm` returns the exact matrix infinity norm; for an ordinary CPU matrix it returns the exact matrix 1-norm. This choice reflects the respective row- and column-oriented storage costs, not equality of the two norms for a general matrix.

A polyadic object is not expanded merely to evaluate a norm. Instead, a block 1-norm estimator probes products with `A` and its adjoint using up to `t` columns and `itmax` iterations. Probe vectors have unit 1-norm, so each resulting column norm is a **lower bound** on the induced matrix 1-norm; the result is the best bound found, not a guaranteed exact norm or upper bound. Sign/adjoint probes seek more informative columns, and the estimate may stop on non-improvement or redundant probes.

## Inputs and outputs

Syntax: `n=cheap_norm(A,t,itmax)`

- `A` — a matrix, or a polyadic representation thereof.
- `t` — optional; positive integer number of probe columns in the polyadic norm estimator; defaults to 1.
- `itmax` — optional; integer maximum number of estimator iterations (at least 2); defaults to 5.
- `n` — infinity-norm for GPU arrays, 1-norm for CPU arrays, and a lower-bound 1-norm estimate for polyadics.

## References

- Higham and Tisseur, algorithm referenced in the source: <https://doi.org/10.1137/S0895479899356080>
- Spinach Wiki page: <https://spindynamics.org/wiki/index.php?title=cheap_norm.m>
