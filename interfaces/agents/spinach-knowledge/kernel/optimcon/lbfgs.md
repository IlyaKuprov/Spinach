# kernel/optimcon/lbfgs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/lbfgs.m`
- Signature: `direction=lbfgs(dx_hist,dg_hist,g)`
- Total lines: 107

## Purpose

Calculates an approximation to the Newton-Raphson search direction for maximising a function using past gradients to build a serviceable substitute to a Hessian. The Hes- sian matrix is not formed explicitly. Syntax: direction=lbfgs(dx_hist,dg_hist,g)

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

- Lines 32-33: Check consistency; implemented by `grumble(dx_hist,dg_hist,g)`.
- Lines 35-36: Curv. tolerance; implemented by `curv_rel_tol=0.01`.
- Lines 38-39: Relevant inner products; implemented by `dgdx=sum(dg_hist.*dx_hist,1)`.
- Lines 43-45: Curvature pair validation; implemented by `valid=isfinite(dgdx)&isfinite(dgdg)&isfinite(dxdx)& (dgdg>0)&(dxdx>0)&(dgdx<-curv_rel_tol*sqrt(dgdg.*dxdx))`.
- Lines 47-48: Dropping bad pairs; implemented by `dx_hist=dx_hist(:,valid)`.
- Lines 51-53: Return steepest ascent if all pairs are bad; implemented by `if isempty(dx_hist)`.
- Lines 57-58: Initialise variables; implemented by `N=size(dx_hist,2)`.
- Lines 62-63: History; implemented by `for n=1:N`.
- Lines 69-70: Scaling of initial Hessian (identity matrix); implemented by `p_k=dg_hist(:,1)'*dx_hist(:,1)/sum(dg_hist(:,1).^2)`.
- Lines 72-73: Direction; implemented by `direction=p_k*g`.
- Lines 79-80: Convert to ascent; implemented by `direction=-direction`.

### Control flow inferred from the code

- Line 53: conditional branch on `isempty(dx_hist)`.
- Line 63: `for` loop over `n=1:N`.
- Line 74: `for` loop over `n=N:-1:1`.

### Key state/data transformations

- Lines 36: computes `curv_rel_tol` using `curv_rel_tol=0.01`.
- Lines 39: computes `dgdx` using `dgdx=sum(dg_hist.*dx_hist,1)`.
- Lines 40: computes `dgdg` using `dgdg=sum(dg_hist.*dg_hist,1)`.
- Lines 41: computes `dxdx` using `dxdx=sum(dx_hist.*dx_hist,1)`.
- Lines 44-45: computes `valid` using `valid=isfinite(dgdx)&isfinite(dgdg)&isfinite(dxdx)& (dgdg>0)&(dxdx>0)&(dgdx<-curv_rel_tol*sqrt(dgdg.*dxdx))`.
- Lines 48: computes `dx_hist` using `dx_hist=dx_hist(:,valid)`.
- Lines 49: computes `dg_hist` using `dg_hist=dg_hist(:,valid)`.
- Lines 54: computes `direction` using `direction=g; return`.
- Lines 58: computes `N` using `N=size(dx_hist,2)`.
- Lines 59: computes `alpha` using `alpha=zeros(1,N)`.
- Lines 60: computes `p` using `p=zeros(1,N)`.
- Lines 64: computes `p(n)` using `p(n)=1/(dg_hist(:,n)'*dx_hist(:,n))`.
- Lines 65: computes `alpha(n)` using `alpha(n)=p(n)*dx_hist(:,n)'*g`.
- Lines 66: computes `g` using `g=g-alpha(n)*dg_hist(:,n)`.
- Lines 70: computes `p_k` using `p_k=dg_hist(:,1)'*dx_hist(:,1)/sum(dg_hist(:,1).^2)`.
- Lines 75: computes `b` using `b=p(n)*dg_hist(:,n)'*direction`.

### Local helper functions

- Line 85: `grumble()` — `function grumble(dx_hist,dg_hist,g_new)`.
  - Representative operation: `if (size(dx_hist,2)~=size(dg_hist,2))||(size(g_new,2)~=1)`.
  - Representative operation: `error('all vector inputs must be column vector arrays.')`.

## Parameters / inputs

- dx_hist -history of x increments, a stack
- of column vectors, from the latest
- to the earliest
- dg_hist -history of gradient increments,
- a stack of column vectors, from
- the latest to the earliest
- g -current gradient

## Outputs

- direction -LBFGS approximation to the
- maximisation step vector

## Implementation structure

- Calculates an approximation to the Newton-Raphson search
- direction for maximising a function using past gradients
- to build a serviceable substitute to a Hessian. The Hes-
- sian matrix is not formed explicitly. Syntax:
- direction=lbfgs(dx_hist,dg_hist,g)
- dx_hist -history of x increments, a stack
- of column vectors, from the latest
- to the earliest
- dg_hist -history of gradient increments,
- a stack of column vectors, from
- the latest to the earliest
- g -current gradient

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `dx_hist()`, `dg_hist()`, `alpha()`, `any()`, `g_new()`.
