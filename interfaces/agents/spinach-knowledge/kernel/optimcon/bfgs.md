# kernel/optimcon/bfgs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/bfgs.m`
- Signature: `H=bfgs(dx_hist,dg_hist,g)`
- Total lines: 135

## Purpose

Calculates a BFGS approximation to the Newton-Raphson search direction for maximising a function using past gradients to build a serviceable substitute to a Hessian. Unlike LBFGS, the pseudo-Hessian matrix is formed explicitly. Syntax: H=bfgs(dx_hist,dg_hist,g)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(dx_hist,dg_hist,g)`.
- Lines 37-38: Curvature tolerance; implemented by `curv_rel_tol=0.01`.
- Lines 40-41: Relevant inner products; implemented by `dgdx=sum(dg_hist.*dx_hist,1)`.
- Lines 45-47: Curvature pair validation; implemented by `valid=isfinite(dgdx)&isfinite(dgdg)&isfinite(dxdx)& (dgdg>0)&(dxdx>0)&(dgdx<-curv_rel_tol*sqrt(dgdg.*dxdx))`.
- Lines 49-50: Dropping bad pairs; implemented by `dx_hist=dx_hist(:,valid)`.
- Lines 53-55: Return identity if all pairs are bad; implemented by `if isempty(dx_hist)`.
- Lines 59-60: Dimension and history length; implemented by `n=numel(g); N=size(dx_hist,2)`.
- Lines 62-63: Pull the first curvature pair; implemented by `dx=dx_hist(:,1); dg=-dg_hist(:,1)`.
- Lines 65-66: Inner products; implemented by `dgdx=dg'*dx; dgdg=dg'*dg`.
- Lines 68-71: Unit matrix scaling; implemented by `if (~isfinite(dgdx))|| (~isfinite(dgdg))|| (dgdx<=eps)||(dgdg<=0)`.
- Lines 73-74: Bad pair; implemented by `gamma=1`.
- Lines 78-79: Good pair; implemented by `gamma=dgdg/dgdx`.
- Lines 83-84: Initial guess; implemented by `H=gamma*eye(n)`.
- Lines 86-87: BFGS updates; implemented by `for k=N:-1:1`.
- Lines 89-90: Pull next curvature pair; implemented by `dx=dx_hist(:,k); dg=-dg_hist(:,k)`.
- Lines 92-93: Make BFGS components; implemented by `dgdx=dg'*dx; Hdx=H*dx; dxHdx=dx'*Hdx`.
- Lines 95-96: Run pair validity checks to protect numerics; implemented by `if (~isfinite(dgdx))||(dgdx<=eps), continue; end`.
- Lines 99-100: Plain BFGS Hessian update; implemented by `H=H-(Hdx*Hdx')/dxHdx+(dg*dg')/dgdx`.

### Control flow inferred from the code

- Line 55: conditional branch on `isempty(dx_hist)`.
- Line 69: conditional branch on `(~isfinite(dgdx))||`.
- Line 87: `for` loop over `k=N:-1:1`.
- Line 96: conditional branch on `(~isfinite(dgdx))||(dgdx<=eps), continue; end`.
- Line 97: conditional branch on `(~isfinite(dxHdx))||(dxHdx<=eps), continue; end`.

### Key state/data transformations

- Lines 38: computes `curv_rel_tol` using `curv_rel_tol=0.01`.
- Lines 41: computes `dgdx` using `dgdx=sum(dg_hist.*dx_hist,1)`.
- Lines 42: computes `dgdg` using `dgdg=sum(dg_hist.*dg_hist,1)`.
- Lines 43: computes `dxdx` using `dxdx=sum(dx_hist.*dx_hist,1)`.
- Lines 46-47: computes `valid` using `valid=isfinite(dgdx)&isfinite(dgdg)&isfinite(dxdx)& (dgdg>0)&(dxdx>0)&(dgdx<-curv_rel_tol*sqrt(dgdg.*dxdx))`.
- Lines 50: computes `dx_hist` using `dx_hist=dx_hist(:,valid)`.
- Lines 51: computes `dg_hist` using `dg_hist=dg_hist(:,valid)`.
- Lines 56: computes `H` using `H=eye(numel(g)); return`.
- Lines 60: computes `n` using `n=numel(g); N=size(dx_hist,2)`.
- Lines 63: computes `dx` using `dx=dx_hist(:,1); dg=-dg_hist(:,1)`.
- Lines 74: computes `gamma` using `gamma=1`.

### Local helper functions

- Line 110: `grumble()` — `function grumble(dx_hist,dg_hist,g_new)`.
  - Representative operation: `if (size(dx_hist,2)~=size(dg_hist,2))||(size(g_new,2)~=1)`.
  - Representative operation: `error('all vector inputs must be column vector arrays.')`.

## Parameters / inputs

- dx_hist -history of x increments, a stack
- of column vectors, from the latest
- to the earliest
- dg_hist -history of gradient increments,
- a stack of column vectors, from
- the latest to the earliest
- g -current gradient (used for sizing)

## Outputs

- H -BFGS approximation to the Hessian
- matrix corresponding to the *nega-
- tive* Hessian of the objective.
- The corresponding ascent directi-
- on is obtained as: direction=H\g

## Implementation structure

- Calculates a BFGS approximation to the Newton-Raphson search
- direction for maximising a function using past gradients to
- build a serviceable substitute to a Hessian. Unlike LBFGS,
- the pseudo-Hessian matrix is formed explicitly. Syntax:
- H=bfgs(dx_hist,dg_hist,g)
- dx_hist -history of x increments, a stack
- of column vectors, from the latest
- to the earliest
- dg_hist -history of gradient increments,
- a stack of column vectors, from
- the latest to the earliest
- g -current gradient (used for sizing)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `dx_hist()`, `dg_hist()`, `any()`, `g_new()`.
