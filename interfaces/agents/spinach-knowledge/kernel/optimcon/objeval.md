# kernel/optimcon/objeval.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/objeval.m`
- Signature: `[data,fx,grad,hess]=objeval(x,objfun_handle,data,spin_system)`
- Total lines: 122

## Purpose

Calls and collect the correct amount of outputs from an objective function -used by optimisation routines. Syntax: [data,fx,grad,hess]=objeval(x,objfun_handle,data,spin_system)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 37-38: Check consistency; implemented by `grumble(objfun_handle,x)`.
- Lines 40-41: Switch between function/gradient/Hessian calls; implemented by `switch nargout`.
- Lines 45-46: Run the objective function for the fidelity; implemented by `[traj_data,fidelity]=feval(objfun_handle,reshape(x,data.x_shape),spin_system)`.
- Lines 48-49: Assign the data structure; implemented by `data.fx_sep_pen=fidelity`.
- Lines 52-53: Compute total error functional; implemented by `fx=fidelity(1)-sum(fidelity(2:end))`.
- Lines 55-56: Increment counters; implemented by `data.count.fx=data.count.fx+1`.
- Lines 60-61: Run the objective function for the fidelity and gradient; implemented by `[traj_data,fidelity,grad]=feval(objfun_handle,reshape(x,data.x_shape),spin_system)`.
- Lines 70-71: Compute the total gradient; implemented by `grad=grad(:,:,1)-sum(grad(:,:,2:end),3); grad=grad(:)`.
- Lines 79-80: Run the objective function for the fidelity, gradient, and Hessian; implemented by `[traj_data,fidelity,grad,hess]=feval(objfun_handle,reshape(x,data.x_shape),spin_system)`.
- Lines 92-93: Compute the total Hessian; implemented by `hess=hess(:,:,1)-sum(hess(:,:,2:end),3)`.
- Lines 102-103: Complain and bomb out; implemented by `error('must have at least two output arguments.')`.

### Control flow inferred from the code

- Line 41: dispatches on `nargout`; cases `2`, `3`, `4`.

### Key state/data transformations

- Lines 46: computes `[traj_data,fidelity]` using `[traj_data,fidelity]=feval(objfun_handle,reshape(x,data.x_shape),spin_system)`.
- Lines 49: computes `data.fx_sep_pen` using `data.fx_sep_pen=fidelity`.
- Lines 50: computes `data.traj_data` using `data.traj_data=traj_data`.
- Lines 53: computes `fx` using `fx=fidelity(1)-sum(fidelity(2:end))`.
- Lines 56: computes `data.count.fx` using `data.count.fx=data.count.fx+1`.
- Lines 61: computes `[traj_data,fidelity,grad]` using `[traj_data,fidelity,grad]=feval(objfun_handle,reshape(x,data.x_shape),spin_system)`.
- Lines 71: computes `grad` using `grad=grad(:,:,1)-sum(grad(:,:,2:end),3); grad=grad(:)`.
- Lines 75: computes `data.count.gfx` using `data.count.gfx=data.count.gfx+1`.
- Lines 80: computes `[traj_data,fidelity,grad,hess]` using `[traj_data,fidelity,grad,hess]=feval(objfun_handle,reshape(x,data.x_shape),spin_system)`.
- Lines 93: computes `hess` using `hess=hess(:,:,1)-sum(hess(:,:,2:end),3)`.
- Lines 98: computes `data.count.hfx` using `data.count.hfx=data.count.hfx+1`.

### Local helper functions

- Line 110: `grumble()` — `function grumble(objfun_handle,x)`. If you can't write it down in English, you can't code it. Peter Halpern
  - Representative operation: `if ~isa(objfun_handle,'function_handle')`.
  - Representative operation: `error('objfun_handle must be a function handle.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `feval()`, `fidelity()`, `grad()`, `hess()`.
