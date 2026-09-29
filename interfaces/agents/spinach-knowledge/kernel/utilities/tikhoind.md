# kernel/utilities/tikhoind.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/tikhoind.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/tikhoind.m)

## Purpose

`tikhoind` computes the analytical Tikhonov-regularised solution to `K*x=y` without any constraints, producing a sign-indefinite output. It minimises `norm(K*x-y,2)^2 + lambda*norm(D*x,2)^2`.

## Behaviour

- Validates inputs via an internal `grumble` consistency check:
  - All inputs (`K`, `D`, `y`, `lam`) must be numeric.
  - The number of rows of `K` must match the number of rows of `y`.
  - `lam` must be a positive real scalar (the check rejects non-real, non-scalar, or negative values).
- Computes the analytical solution as `x = ((K'*K) + lam*(D'*D)) \ (K'*y)`.
- If more than one output is requested, computes `err = norm(K*x-y,2)^2`.
- If more than two outputs are requested, computes `reg = norm(D*x,2)^2`.
- The kernel matrix `K` may be complex and non-square; `y` may be complex.
- For best numerical performance, the source recommends scaling `K` to have approximately unit 2-norm and `y` to have approximately unit 1-norm.
- The source notes that `tikhonov.m` provides the positive-constrained solver.

## Inputs and outputs

**Syntax:** `[x,err,reg]=tikhoind(K,D,y,lam)`

**Inputs:**

- `K` — kernel matrix, may be complex, may be non-square.
- `D` — regularisation matrix.
- `y` — column vector, may be complex.
- `lam` — Tikhonov regularisation parameter; must be a positive real scalar.

**Outputs:**

- `x` — a real vector, a minimum of `norm(K*x-y,2)^2 + lambda*norm(D*x,2)^2`.
- `err` — error signal, `norm(K*x-y,2)^2`.
- `reg` — regularisation signal, `norm(D*x,2)^2`.

## References

- Spinach Wiki: [tikhoind.m](https://spindynamics.org/wiki/index.php?title=tikhoind.m)
