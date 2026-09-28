# kernel/optimcon/bracketing.m

- Signature: `[a,b,alpha,fx,gfx,next_act,data]=bracketing(cost_function,alpha,dir,x_0,...`

## Purpose

Expands a trial step to find a bracket containing an acceptable line-search point, or accepts the step if the Wolfe tests are met before sectioning is needed.

## Syntax

```matlab
[a,b,alpha,fx,gfx,next_act,data]=bracketing(cost_function,alpha,dir,x_0,...
                                              fx_0,gfx_0,data,spin_system)
```

## Inputs

- `cost_function` — objective-function handle.
- `alpha` — initial trial step length.
- `dir` — search direction.
- `x_0` — current optimisation vector.
- `fx_0` — objective value at `x_0`.
- `gfx_0` — gradient at `x_0`.
- `data` — optimisation workspace passed to and updated by the objective evaluation.
- `spin_system` — Spinach data structure containing line-search settings and any frozen-coordinate mask.

## Outputs

- `a`, `b` — lower and upper bracket-point structures.
- `alpha` — current step length; it is the accepted length when the Wolfe conditions are met.
- `fx`, `gfx` — objective value and gradient at the current trial point.
- `next_act` — continuation tag: `'sectioning'` when sectioning should follow, or `'none'` when no further line-search stage is needed.
- `data` — updated optimisation workspace.

## Algorithm

Frozen coordinates are removed from the search direction and initial gradient when a freeze mask is present. The routine evaluates trial points while expanding or tightening the bracket, and uses cubic interpolation within the bracket to select subsequent trial steps.
