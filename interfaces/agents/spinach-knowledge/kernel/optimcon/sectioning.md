# kernel/optimcon/sectioning.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/sectioning.m`
- Signature: `[alpha,fx_1,gfx_1,exitflag,data]=sectioning(cost_function,a,b,x_0,fx_0,...`
- Total lines: 182

## Purpose

Refines a previously found step bracket by repeated cubic interpolation until a step satisfying Wolfe tests is found or the bracket collapses below numerical accuracy. Syntax: [alpha,fx_1,gfx_1,exitflag,data]=... sectioning(cost_function,A,B,x_0,fx_0,gfx_0,... dir,data,spin_system)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 52-53: Check consistency; implemented by `grumble(cost_function,a,b,x_0,fx_0,gfx_0,dir)`.
- Lines 55-56: Apply coordinate freezing mask when requested; implemented by `if ~isempty(spin_system.control.freeze)`.
- Lines 61-62: Iterate until Wolfe conditions pass or bracket collapses; implemented by `while true`.
- Lines 64-65: Build reduced interpolation bounds inside the bracket; implemented by `end_a=a.alpha+spin_system.control.ls_tau2*(b.alpha-a.alpha)`.
- Lines 68-70: Maximise cubic model inside reduced interpolation bounds; implemented by `alpha=cubic_interp(end_a,end_b,a.alpha,b.alpha, a.fx,a.gfx'*dir,b.fx,b.gfx'*dir)`.
- Lines 72-73: Stop when interpolation displacement is numerically unresolved; implemented by `if abs((alpha-a.alpha)*(a.gfx'*dir))<eps(max(1,abs(a.fx)))`.
- Lines 78-79: Evaluate objective and gradient at current trial step; implemented by `[data,fx_1,gfx_1]=objeval(x_0+alpha*dir,cost_function,data,spin_system)`.
- Lines 81-82: Apply coordinate freezing mask to the new gradient; implemented by `if ~isempty(spin_system.control.freeze)`.
- Lines 86-87: Store current lower bracket endpoint before updates; implemented by `tmp=a`.
- Lines 89-91: Update bracket based on Armijo and monotonicity tests; implemented by `if (~alpha_conds(1,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system))|| (~alpha_conds(0,alpha,a.fx,fx_1,a.gfx,gfx_1,dir,spin_system))`.
- Lines 93-94: Move upper bracket endpoint to current trial step; implemented by `b.alpha=alpha; b.fx=fx_1; b.gfx=gfx_1`.
- Lines 98-99: Accept trial step if curvature condition is satisfied; implemented by `if alpha_conds(2,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system)`.
- Lines 103-104: Move lower bracket endpoint to current trial step; implemented by `a.alpha=alpha; a.fx=fx_1; a.gfx=gfx_1`.
- Lines 106-107: Swap upper endpoint when derivative sign condition fails; implemented by `if (a.alpha-b.alpha)*(gfx_1'*dir)>=0`.
- Lines 113-114: Terminate when bracket width is numerically unresolved; implemented by `if abs((b.alpha-a.alpha)*(a.gfx'*dir))<eps(max(1,abs(a.fx)))`.

### Control flow inferred from the code

- Line 56: conditional branch on `~isempty(spin_system.control.freeze)`.
- Line 62: `while` loop over `true`.
- Line 73: conditional branch on `abs((alpha-a.alpha)*(a.gfx'*dir))<eps(max(1,abs(a.fx)))`.
- Line 82: conditional branch on `~isempty(spin_system.control.freeze)`.
- Line 90: conditional branch on `(~alpha_conds(1,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system))||`.
- Line 99: conditional branch on `alpha_conds(2,alpha,fx_0,fx_1,gfx_0,gfx_1,dir,spin_system)`.
- Line 107: conditional branch on `(a.alpha-b.alpha)*(gfx_1'*dir)>=0`.
- Line 114: conditional branch on `abs((b.alpha-a.alpha)*(a.gfx'*dir))<eps(max(1,abs(a.fx)))`.

### Key state/data transformations

- Lines 57: computes `dir` using `dir=dir.*(~spin_system.control.freeze(:))`.
- Lines 58: computes `gfx_0` using `gfx_0=gfx_0.*(~spin_system.control.freeze(:))`.
- Lines 69-70: computes `alpha` using `alpha=cubic_interp(end_a,end_b,a.alpha,b.alpha, a.fx,a.gfx'*dir,b.fx,b.gfx'*dir)`.
- Lines 74: computes `exitflag` using `exitflag=-2; alpha=a.alpha`.
- Lines 75: computes `fx_1` using `fx_1=a.fx; gfx_1=a.gfx; return`.
- Lines 79: computes `[data,fx_1,gfx_1]` using `[data,fx_1,gfx_1]=objeval(x_0+alpha*dir,cost_function,data,spin_system)`.
- Lines 83: computes `gfx_1` using `gfx_1=gfx_1.*(~spin_system.control.freeze(:))`.
- Lines 87: computes `tmp` using `tmp=a`.
- Lines 94: computes `b.alpha` using `b.alpha=alpha; b.fx=fx_1; b.gfx=gfx_1`.
- Lines 104: computes `a.alpha` using `a.alpha=alpha; a.fx=fx_1; a.gfx=gfx_1`.
- Lines 108: computes `b` using `b=tmp`.

### Local helper functions

- Line 124: `grumble()` — `function grumble(cost_function,a,b,x_0,fx_0,gfx_0,dir)`.
  - Representative operation: `if ~isa(cost_function,'function_handle')`.
  - Representative operation: `error('cost_function must be a function handle.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `cubic_interp()`, `eps()`, `objeval()`, `alpha_conds()`, `isstruct()`, `isfield()`, `isscalar()`, `iscolumn()`, `isequal()`.
