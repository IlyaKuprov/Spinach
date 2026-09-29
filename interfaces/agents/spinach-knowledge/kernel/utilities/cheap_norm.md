# kernel/utilities/cheap_norm.m

## Purpose

Returns the cheapest available matrix norm depending on how the matrix is represented: the infinity-norm for GPU arrays, the 1-norm for CPU arrays, and a lower-bound 1-norm estimate for polyadic objects, which can only multiply vectors. The header notes that some norms are vastly more expensive than others and that this function uses the cheapest ones available.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/cheap_norm.m>

## Behaviour

- Defaults are applied when arguments are omitted: `t` defaults to 1 and `itmax` defaults to 5.
- Input validation (`grumble`) requires `A` to be numeric (matrix or polyadic), `t` a positive real integer, and `itmax` a real integer greater than one; violations raise errors.
- If `A` is a `gpuArray`, the function returns `norm(A,inf)` immediately.
- If `A` is not a `polyadic` object (CPU case), the function returns `norm(A,1)` immediately.
- For polyadic `A`, the function runs a block 1-norm estimator based on Algorithm 2.4 of Higham and Tisseur's paper:
  - `t` is capped at the column dimension of `A` (`t=min(t,col_dim)`).
  - The initial probe matrix `X` has `t` columns: the first is all ones, the rest are random ±1 vectors; columns are resampled (up to 100 tries each) until no pair of columns is parallel, then `X` is divided by the column dimension.
  - Each iteration computes `Y=A*X`, takes column sums of absolute values as candidate estimates, and tracks the best estimate and its column index.
  - The loop terminates early when the estimate does not improve, when `itmax` iterations are exceeded, when all sign probes are redundant (real case), when the best column has been reached via adjoint products `Z=A'*S`, or when all candidate unit-vector probes have already been used.
  - The phase matrix `S=sign(Y)` uses the paper's zero convention, replacing zeros with 1.
  - In the real case, sign columns are resampled (up to 100 tries) to avoid parallelism with other sign columns and the previous phase matrix; exceeding 100 tries returns the current best estimate.
  - Subsequent probe matrices are built from unit vectors at row indices with the largest scores from `max(abs(Z),[],2)`, sorted in descending order, with previously used indices deprioritised.
  - The function returns the best lower bound encountered (`est_old`).

## Inputs and outputs

Syntax: `n=cheap_norm(A,t,itmax)`

- `A` — a matrix, or a polyadic representation thereof.
- `t` — optional; number of probe columns in the polyadic norm estimator; defaults to 1.
- `itmax` — optional; maximum number of estimator iterations; defaults to 5.
- `n` — infinity-norm for GPU arrays, 1-norm for CPU arrays, and a lower-bound 1-norm estimate for polyadics.

## References

- Higham and Tisseur, algorithm referenced in the source: <https://doi.org/10.1137/S0895479899356080>
- Spinach Wiki page: <https://spindynamics.org/wiki/index.php?title=cheap_norm.m>
