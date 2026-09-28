# kernel/grids/gaussleg.m

- Signature: `[x,w]=gaussleg(a,b,n)`

## Purpose

Computes Gauss-Legendre points and weights on the interval `[a,b]`. The requested accuracy order `n` produces `n+1` points.

## Parameters / inputs

- `a` — finite real scalar left endpoint, with `a < b`.
- `b` — finite real scalar right endpoint.
- `n` — positive real integer accuracy order; values above 40 are rejected.

## Outputs

- `x` — Gauss-Legendre points on `[a,b]`, sorted in ascending order.
- `w` — corresponding Gauss-Legendre weights, reordered with `x`.

## Method

The routine initializes nodes in `[-1,1]`, refines them by Newton-Raphson iteration, computes the weights, maps the nodes to `[a,b]`, and sorts the points and their weights together.

Source reference: <https://spindynamics.org/wiki/index.php?title=gaussleg.m>
