# kernel/optimcon/bfgs_upd.m

- Signature: `H=bfgs_upd(H,dx,dg)`

## Purpose

Performs one dense BFGS update for maximisation. `H` approximates the negative Hessian of the objective, and `dg` is the gradient increment between steps.

## Algorithm

The update uses the sign-adjusted gradient increment `-dg`. A curvature safeguard rejects pairs that are non-finite or do not satisfy the required negative-curvature test on `dg' * dx`. If `H` is empty and the pair is rejected, the routine returns an identity matrix; if a supplied `H` is paired with a rejected step, it is symmetrised and returned unchanged. For an empty `H` with an accepted pair, a scaled identity initializes the approximation before the BFGS update.

## Syntax

```matlab
H=bfgs_upd(H,dx,dg)
```

## Inputs

- `H` — current approximation to the negative Hessian, or `[]` on the first call.
- `dx` — argument increment between the current and previous steps.
- `dg` — gradient increment between the current and previous steps.

## Output

- `H` — updated BFGS approximation to the negative Hessian.
