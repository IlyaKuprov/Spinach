# kernel/optimcon/sectioning.m

- Signature: `[alpha,fx_1,gfx_1,exitflag,data]=sectioning(cost_function,a,b,x_0,fx_0,...`

## Purpose

Refines a previously found step bracket by repeated cubic interpolation until a step satisfying Wolfe tests is found or the bracket collapses below numerical accuracy. Syntax: [alpha,fx_1,gfx_1,exitflag,data]=... sectioning(cost_function,A,B,x_0,fx_0,gfx_0,... dir,data,spin_system)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

## Parameters / inputs

- cost_function -objective function handle
- a -lower bracket structure with
- fields alpha, fx, and gfx
- b -upper bracket structure with
- fields alpha, fx, and gfx
- x_0 -current optimisation vector
- fx_0 -objective value at x_0
- gfx_0 -gradient at x_0
- dir -search direction vector
- data -optimisation workspace structure
- spin_system -Spinach data structure with
- sectioning settings

## Outputs

- alpha -accepted step length
- fx_1 -objective value at alpha
- gfx_1 -gradient at alpha
- exitflag -0 on success, -2 on failure
- data -updated optimisation workspace

## Implementation structure

- Refines a previously found step bracket by repeated cubic
- interpolation until a step satisfying Wolfe tests is found
- or the bracket collapses below numerical accuracy. Syntax:
- [alpha,fx_1,gfx_1,exitflag,data]=...
- sectioning(cost_function,A,B,x_0,fx_0,gfx_0,...
- dir,data,spin_system)
- cost_function -objective function handle
- a -lower bracket structure with
- fields alpha, fx, and gfx
- b -upper bracket structure with
- x_0 -current optimisation vector
- fx_0 -objective value at x_0
