# kernel/utilities/iseye.m

## Purpose

Returns `true` for unit (identity) matrices. The test is designed to be computationally affordable. Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/iseye.m>

## Behaviour

- Syntax: `verdict=iseye(M)`.
- Consistency is enforced first: if `M` is not numeric, the function errors with `'M must be numeric.'`.
- If `M` is not square, `verdict` is `false`.
- Otherwise, if `M` is not diagonal (`~isdiag(M)`), `verdict` is `false`.
- Otherwise, a random test vector `a=randn(size(M,2),1)` is generated, and `M*a` is compared with `a` via `nnz(M*a-a)`. If `nnz(M*a-a)~=0`, `verdict` is `false`; otherwise `verdict` is `true`.

## Inputs and outputs

**Inputs**

- `M` — a matrix (must be numeric).

**Outputs**

- `verdict` — `true` or `false`.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=iseye.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/iseye.m>
