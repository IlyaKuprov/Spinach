# kernel/optimcon/objeval.m

- Signature: `[data,fx,grad,hess]=objeval(x,objfun_handle,data,spin_system)`

## Purpose

Adapts an objective-function call to the number of outputs requested by the optimisation routine. It reshapes the parameter vector using `data.x_shape`, calls the objective with `spin_system`, and combines fidelity, gradient, and Hessian contributions into the total objective and its derivatives.

## Parameters / inputs

- `x` — non-empty real numeric parameter vector.
- `objfun_handle` — objective function handle; it must provide at least two outputs.
- `data` — optimisation state, including `x_shape` and counters; objective diagnostics and trajectory data are stored here.
- `spin_system` — passed to the objective function.

## Outputs

- `data` — updated state containing the separate fidelity contributions, trajectory data, and evaluation counters.
- `fx` — total objective, computed as the first fidelity contribution minus the sum of the remaining contributions.
- `grad` — combined gradient, returned when three or more outputs are requested.
- `hess` — combined Hessian, returned when four outputs are requested.

The function supports requests for two, three, or four outputs; other output counts raise an error. The source marks the function for elimination in a future release.

[Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=objeval.m)
