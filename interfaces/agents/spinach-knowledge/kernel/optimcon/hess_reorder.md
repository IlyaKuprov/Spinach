# kernel/optimcon/hess_reorder.m

- Signature: `hess=hess_reorder(hess,K,N)`

## Purpose

Reorders a Hessian whose variables are laid out control-channel first and time-point second to the matching time-point-first, control-channel-second ordering. This converts the ordering `[X1 Y1 Z1 X2 Y2 Z2 ...]` to `[X1 X2 ... Y1 Y2 ... Z1 Z2 ...]` for control channels X, Y, Z.

## Parameters / inputs

- `hess` — square `(K*N)	imes(K*N)` Hessian in control-channel-first ordering.
- `K` — positive integer number of control channels.
- `N` — positive integer number of time points.

## Output

- `hess` — the same Hessian reordered to time-point-first ordering.

## Implementation

The function reshapes the matrix to `[K N K N]`, permutes both variable axes with `[2 1 4 3]`, then reshapes it back to `[N*K N*K]`.

[Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=hess_reorder.m)
