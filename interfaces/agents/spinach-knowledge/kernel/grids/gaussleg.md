# kernel/grids/gaussleg.m

- Signature: `[x,w]=gaussleg(a,b,n)`

## Purpose

Constructs a Gauss-Legendre quadrature rule on the finite interval `[a,b]`.

## Rule, dimensions, and reproducibility

The requested positive integer `n` gives `n+1` nodes: the routine refines the Legendre-polynomial roots on `[-1,1]`, computes their Gauss-Legendre weights, maps both to `[a,b]`, then sorts the nodes in ascending order with their corresponding weights. Thus `x` and `w` are matching `(n+1)`-element column vectors. As a Gauss-Legendre rule with `n+1` nodes, it integrates polynomials through degree `2*n+1` exactly in exact arithmetic.

Node starts are computed from a fixed formula and refined by Newton iteration; no random sampling is used. Iteration stops when the largest node update is no greater than machine epsilon. The implementation rejects `n>40` and advises subdividing the interval instead.

## Parameters / inputs

- `a`, `b` — finite real scalar endpoints with `a < b`.
- `n` — finite positive real integer no greater than 40.

## Outputs

- `x` — ascending Gauss-Legendre nodes on `[a,b]`.
- `w` — corresponding integration weights for `[a,b]`.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/gaussleg.m)
<https://spindynamics.org/wiki/index.php?title=gaussleg.m>
