# kernel/optimcon/bfgs_upd.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/bfgs_upd.m`
- Signature: `H=bfgs_upd(H,dx,dg)`
- Total lines: 118

## Purpose

Performs a one-step BFGS Hessian update for maximisation using the argument and gradient increments from the previous step.

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- This routine performs a single dense BFGS Hessian update for a maximisation problem. The sign conventions matter: the stored matrix approximates the negative Hessian so that solving H\g yields an ascent direction.
- The curvature test rejects bad secant pairs when dg^T dx does not have the sign and magnitude expected for locally concave behaviour. That protects the update from producing indefinite or numerically meaningless curvature models.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.
- If no valid curvature information exists yet, the code scales an identity matrix using y^T y / y^T dx, a standard quasi-Newton initialisation that roughly matches curvature along the first accepted step.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(H,dx,dg)`.
- Lines 35-36: Increments into columns; implemented by `dx=dx(:); dg=dg(:); y=-dg`.
- Lines 38-39: Curvature tolerance; implemented by `curv_rel_tol=0.01`.
- Lines 41-42: Relevant inner products; implemented by `dgdx=dg'*dx; dgdg=dg'*dg; dxdx=dx'*dx`.
- Lines 44-46: Curvature pair validation; implemented by `pair_ok=isfinite(dgdx)&&isfinite(dgdg)&&isfinite(dxdx)&& (dgdg>0)&&(dxdx>0)&&(dgdx<-curv_rel_tol*sqrt(dgdg*dxdx))`.
- Lines 48-49: Return identity matrix if first pair is bad; implemented by `if isempty(H)&&(~pair_ok), H=eye(numel(dx)); return; end`.
- Lines 51-52: Return unchanged if pair is bad; implemented by `if (~isempty(H))&&(~pair_ok), H=real((H+H')/2); return; end`.
- Lines 54-55: Initialise Hessian estimate when absent; implemented by `if isempty(H)`.
- Lines 57-58: Relevant inner products; implemented by `ydx=y'*dx; yyy=y'*y`.
- Lines 60-63: Unit matrix scaling; implemented by `if (~isfinite(ydx))|| (~isfinite(yyy))|| (ydx<=eps)||(yyy<=0)`.
- Lines 69-70: Initial Hessian guess; implemented by `H=gamma*eye(numel(dx))`.
- Lines 74-75: Build BFGS components; implemented by `ydx=y'*dx; Hdx=H*dx; dxHdx=dx'*Hdx`.
- Lines 77-78: Return unchanged if unsafe; implemented by `if (~isfinite(ydx))||(ydx<=eps)`.
- Lines 85-86: Plain BFGS Hessian update; implemented by `H=H-(Hdx*Hdx')/dxHdx+(y*y')/ydx`.
- Lines 88-89: Numerical clean-up; implemented by `H=real((H+H')/2)`.

### Control flow inferred from the code

- Line 49: conditional branch on `isempty(H)&&(~pair_ok), H=eye(numel(dx)); return; end`.
- Line 52: conditional branch on `(~isempty(H))&&(~pair_ok), H=real((H+H')/2); return; end`.
- Line 55: conditional branch on `isempty(H)`.
- Line 61: conditional branch on `(~isfinite(ydx))||`.
- Line 78: conditional branch on `(~isfinite(ydx))||(ydx<=eps)`.
- Line 81: conditional branch on `(~isfinite(dxHdx))||(dxHdx<=eps)`.

### Key state/data transformations

- Lines 36: computes `dx` using `dx=dx(:); dg=dg(:); y=-dg`.
- Lines 39: computes `curv_rel_tol` using `curv_rel_tol=0.01`.
- Lines 42: computes `dgdx` using `dgdx=dg'*dx; dgdg=dg'*dg; dxdx=dx'*dx`.
- Lines 45-46: computes `pair_ok` using `pair_ok=isfinite(dgdx)&&isfinite(dgdg)&&isfinite(dxdx)&& (dgdg>0)&&(dxdx>0)&&(dgdx<-curv_rel_tol*sqrt(dgdg*dxdx))`.
- Lines 58: computes `ydx` using `ydx=y'*dx; yyy=y'*y`.
- Lines 64: computes `gamma` using `gamma=1`.
- Lines 70: computes `H` using `H=gamma*eye(numel(dx))`.

### Local helper functions

- Line 94: `grumble()` — `function grumble(H,dx,dg)`.
  - Representative operation: `if (~isnumeric(dx))||(~isreal(dx))||(~isvector(dx))||isempty(dx)`.
  - Representative operation: `error('dx must be a non-empty real numeric vector.')`.

## Syntax

```matlab
H=bfgs_upd(H,dx,dg)
```

## Parameters / inputs

- H -current BFGS approximation to the Hessian
- matrix corresponding to the *negative*
- Hessian of the objective; use [] on the
- first call
- dx -increment in arguments between the current
- and the previous step
- dg -increment in gradients between the current
- and the previous step

## Outputs

- H -updated BFGS approximation to the Hessian
- matrix corresponding to the *negative*
- Hessian of the objective

## Implementation structure

- Performs a one-step BFGS Hessian update for maximisation using
- the argument and gradient increments from the previous step.
- H=bfgs_upd(H,dx,dg)
- H -current BFGS approximation to the Hessian
- matrix corresponding to the *negative*
- Hessian of the objective; use [] on the
- first call
- dx -increment in arguments between the current
- and the previous step
- dg -increment in gradients between the current
- H -updated BFGS approximation to the Hessian
- Hessian of the objective

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isvector()`.
