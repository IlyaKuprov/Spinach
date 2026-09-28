# kernel/optimcon/bfgs.m

- Signature: `H=bfgs(dx_hist,dg_hist,g)`

## Purpose

Constructs a dense BFGS approximation to the negative Hessian of an objective being maximised. The resulting ascent direction is obtained by solving `H\g`.

## Algorithm

History columns are considered from newest to oldest. The routine rejects non-finite or insufficient-curvature step/gradient pairs, initializes the matrix scale from the first retained pair, and applies the BFGS updates to the remaining valid pairs. If no pair survives the curvature checks, it returns the identity matrix. This is a full-matrix BFGS approximation, not a limited-memory update.

## Syntax

```matlab
H=bfgs(dx_hist,dg_hist,g)
```

## Inputs

- `dx_hist` — history of argument increments, with one column per pair, newest first.
- `dg_hist` — corresponding gradient increments, in the same column order.
- `g` — current gradient column vector; its length determines the matrix dimension.

## Output

- `H` — dense BFGS approximation to the negative objective Hessian; use `H\g` for the corresponding Newton-like ascent direction.
