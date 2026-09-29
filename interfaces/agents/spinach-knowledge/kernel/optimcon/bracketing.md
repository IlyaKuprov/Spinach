# kernel/optimcon/bracketing.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/optimcon/bracketing.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=bracketing.m)

## Purpose

For an objective being maximised, evaluate a trial step and either accept it under the line-search tests or return a bracket for a later sectioning stage. This is a line-search helper, not a general constrained optimiser.

## Syntax

`[a,b,alpha,fx,gfx,next_act,data]=bracketing(cost_function,alpha,dir,x_0,fx_0,gfx_0,data,spin_system)`

## Inputs

- `cost_function` — objective function handle evaluated by `objeval`.
- `alpha` — scalar trial step multiplier. The input check requires a nonempty numeric real scalar; it does not impose a positivity or finiteness check.
- `dir` — nonempty real numeric search-direction column vector.
- `x_0` — nonempty real numeric starting-point column vector.
- `fx_0` — nonempty real numeric scalar objective value at `x_0`.
- `gfx_0` — nonempty real numeric gradient column vector at `x_0`.
- `data` — workspace passed through `objeval` and returned updated; this routine does not validate it.
- `spin_system` — required state used for control settings. The implementation reads `spin_system.control.freeze` and `spin_system.control.ls_tau1`, but does not validate the state structure.

`dir`, `x_0`, and `gfx_0` must have matching dimensions. These vectors and `fx_0` are checked for numeric, real, and shape requirements but not for finiteness. The scalar `alpha` is not checked for positivity.

## Outputs

- `a` and `b` — endpoint records with fields `alpha`, `fx`, and `gfx`; those fields are populated when a bracket is captured and otherwise may be empty.
- `alpha`, `fx`, and `gfx` — current trial step, objective, and gradient.
- `next_act` — `'sectioning'` when a bracket is handed to sectioning, or `'none'` when the line-search step is accepted without that stage.
- `data` — updated evaluator workspace.

## Search behaviour

The first trial point is `x_0+alpha*dir`. The routine evaluates its objective and gradient, then applies Armijo, monotonicity, curvature, and directional-derivative-sign tests through `alpha_conds`. A failed sufficient-increase or monotonicity test, or a change in directional-derivative sign, records bracket endpoints and returns `next_act='sectioning'`. If the curvature test passes, it returns the accepted point with `next_act='none'`. Otherwise it advances the trial using cubic interpolation over a window extended using `spin_system.control.ls_tau1`. Non-finite expansion bounds, objective values, or directional derivatives raise an error describing apparent unbounded increase.

If `spin_system.control.freeze` is nonempty, frozen coordinates are removed from the search direction and masked out of the initial and trial gradients. No other constraint handling is performed here. `alpha` is a scalar multiplier; the source specifies no universal units for the optimisation coordinates, objective, or gradient.
