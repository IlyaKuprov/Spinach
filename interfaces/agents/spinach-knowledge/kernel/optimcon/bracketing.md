# kernel/optimcon/bracketing.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/bracketing.m`
- Signature: `[a,b,alpha,fx,gfx,next_act,data]=bracketing(cost_function,alpha,dir,x_0,...`
- Total lines: 177

## Purpose

Expands a trial step into a bracket that contains an acceptable line search point, or accepts the step directly if the Wolfe tests are met before sectioning becomes necessary.

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Syntax

```matlab
[A,B,alpha,fx,gfx,next_act,data]=...
bracketing(cost_function,alpha,dir,x_0,fx_0,...
gfx_0,data,spin_system)
```

## Parameters / inputs

- cost_function -objective function handle
- alpha -initial trial step length
- dir -search direction vector
- x_0 -current optimisation vector
- fx_0 -objective value at x_0
- gfx_0 -gradient at x_0
- data -optimisation workspace structure
- spin_system -Spinach data structure with
- line search settings

## Outputs

- A -lower bracket point structure
- B -upper bracket point structure
- alpha -accepted step length when found
- fx -objective value at accepted step
- gfx -gradient at accepted step
- next_act -continuation tag, either
- 'sectioning' or 'none'
- data -updated optimisation workspace

## Implementation structure

- Expands a trial step into a bracket that contains an acceptable
- line search point, or accepts the step directly if the Wolfe
- tests are met before sectioning becomes necessary.
- [A,B,alpha,fx,gfx,next_act,data]=...
- bracketing(cost_function,alpha,dir,x_0,fx_0,...
- gfx_0,data,spin_system)
- cost_function -objective function handle
- alpha -initial trial step length
- dir -search direction vector
- x_0 -current optimisation vector
- fx_0 -objective value at x_0
- gfx_0 -gradient at x_0

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `objeval()`, `alpha_conds()`, `all()`, `cubic_interp()`, `isscalar()`, `iscolumn()`, `isequal()`.
