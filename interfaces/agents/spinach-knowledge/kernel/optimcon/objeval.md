# kernel/optimcon/objeval.m

- Signature: `[data,fx,grad,hess]=objeval(x,objfun_handle,data,spin_system)`

## Purpose

Calls and collect the correct amount of outputs from an objective function -used by optimisation routines. Syntax: [data,fx,grad,hess]=objeval(x,objfun_handle,data,spin_system)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

## Parameters / inputs

- x -objective function argument
- objfun_handle -handle to the objective function
- data -data structure inherited from
- fmaxnewton.m

## Outputs

- data -modified data structure with
- diagnostics from the objective
- fx -objective function value at x
- grad -gradient of the objective function
- at x
- hess -Hessian of the objective function
- at x
- Note: this function will be eliminated in a future release.

## Implementation structure

- Calls and collect the correct amount of outputs from an objective
- function -used by optimisation routines. Syntax:
- [data,fx,grad,hess]=objeval(x,objfun_handle,data,spin_system)
- x -objective function argument
- objfun_handle -handle to the objective function
- data -data structure inherited from
- fmaxnewton.m
- data -modified data structure with
- diagnostics from the objective
- fx -objective function value at x
- grad -gradient of the objective function
- at x
