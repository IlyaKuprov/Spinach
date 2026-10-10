# kernel/optimcon/bfgs_upd.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/bfgs_upd.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=bfgs_upd.m)

## Purpose

Apply one dense BFGS update for maximisation. `H` approximates the negative objective Hessian; `dx` and `dg` are argument and gradient increments. The update uses the sign-adjusted gradient increment `y=-dg`. It updates curvature information only; it does not evaluate the objective or impose constraints.

## Syntax

`H=bfgs_upd(H,dx,dg)`

## Inputs

- `H` — existing real square approximation, or `[]` for initialisation.
- `dx` — nonempty real numeric vector of argument increments.
- `dg` — nonempty real numeric vector of gradient increments, with the same number of elements as `dx`.

The implementation accepts row or column vectors and reshapes both increments into columns. It checks that a supplied `H` is real, square, and dimensionally compatible. It does not require the increments or `H` to be finite at input validation; non-finite increment curvature fails the pair test.

## Update and output

- `H` — updated real symmetric approximation to the negative objective Hessian.

A curvature pair is used only when the finite inner products are positive in norm and `dg' * dx < -0.01*norm(dg)*norm(dx)`. With an empty `H`, a rejected pair returns an identity matrix sized to `dx`; a usable pair initialises a scaled identity and is then applied in the same call. With an existing matrix, a rejected pair leaves its symmetrised value unchanged. The BFGS update also returns that symmetrised value without updating if its denominators are non-finite or no larger than machine `eps`.

The routine assigns no physical units; increments and gradients retain the caller's optimisation-coordinate and objective units.
