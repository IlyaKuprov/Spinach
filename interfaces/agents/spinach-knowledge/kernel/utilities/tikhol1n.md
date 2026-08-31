# kernel/utilities/tikhol1n.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/tikhol1n.m`
- Signature: `[x,err,reg]=tikhol1n(A,y,nnzt)`
- Total lines: 170

## Purpose

L1 norm Tikhonov regularised solver for A*x=y where A is an ill-conditioned matrix. The error functional is norm(A*x-y,2)^2+lambda*norm(x,1), it is minimised using the FISTA algorithm. The user specifies the de- sired number of non-zeroes, lambda parameter is then found by bracketing / bisection. Syntax: [x,err,reg]=tikhol1n(A,y,nnzt)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(A,y,nnzt)`.
- Lines 38-39: Tolerances; implemented by `normest_tol=1e-3`.
- Lines 42-43: Pre-compute CT; implemented by `A_ct=ctranspose(A)`.
- Lines 45-46: Lipschitz constant and initial threshold; implemented by `L=2*(1+normest_tol)*normest(A,normest_tol)^2; thr=1/L`.
- Lines 48-49: Complex soft thresholding function; implemented by `soft_thr=@(x,thr)sign(x).*max(abs(x)-thr,0)`.
- Lines 51-52: Zero initial guess; implemented by `x=zeros(size(A,2),1); x_old=x`.
- Lines 54-55: Threshold brackets for ZF targeting; implemented by `thr_lower=0; thr_upper=max(abs(A_ct*y))`.
- Lines 57-58: Iteration counters and momentum; implemented by `t=1; iter_count=0; converged=false`.
- Lines 60-61: Tell the user we are starting to iterate; implemented by `disp(['FISTA called with target nnz = ' int2str(nnzt)])`.
- Lines 63-64: FISTA loop; implemented by `while ~converged`.
- Lines 66-67: Error vector; implemented by `err_vec=A*x-y`.
- Lines 69-70: Error gradient; implemented by `g=2*(A_ct*err_vec)`.
- Lines 72-73: Compute proximal point; implemented by `x_prox=soft_thr(x-g/L,thr)`.
- Lines 75-76: Update Nesterov momentum; implemented by `t_new=0.5*(1+sqrt(1+4*t^2))`.
- Lines 78-79: Take the step; implemented by `x=x_prox+((t-1)/t_new)*(x_prox-x_old)`.
- Lines 81-82: Check convergence; implemented by `soln_norm=norm(x_prox,2)`.
- Lines 86-87: Close the loop; implemented by `t=t_new; x_old=x_prox`.
- Lines 90-91: Progress report and nnz targeting; implemented by `if mod(iter_count,1000)==0 || converged`.

### Control flow inferred from the code

- Line 64: `while` loop over `~converged`.
- Line 91: conditional branch on `mod(iter_count,1000)==0 || converged`.
- Line 94: conditional branch on `nnz(x)==0 || (converged && ~ismember(nnz(x),[nnzt-1, nnzt, nnzt+1]))`.
- Line 97: conditional branch on `nnz(x)>nnzt`.
- Line 120: conditional branch on `(thr_upper-thr_lower)/(thr_upper+thr_lower)<1e-6`.

### Key state/data transformations

- Lines 39: computes `normest_tol` using `normest_tol=1e-3`.
- Lines 40: computes `step_norm_tol` using `step_norm_tol=1e-6`.
- Lines 43: computes `A_ct` using `A_ct=ctranspose(A)`.
- Lines 46: computes `L` using `L=2*(1+normest_tol)*normest(A,normest_tol)^2; thr=1/L`.
- Lines 49: computes `soft_thr` using `soft_thr=@(x,thr)sign(x).*max(abs(x)-thr,0)`.
- Lines 52: computes `x` using `x=zeros(size(A,2),1); x_old=x`.
- Lines 55: computes `thr_lower` using `thr_lower=0; thr_upper=max(abs(A_ct*y))`.
- Lines 58: computes `t` using `t=1; iter_count=0; converged=false`.
- Lines 61: computes `disp(['FISTA called with target nnz` using `disp(['FISTA called with target nnz = ' int2str(nnzt)])`.
- Lines 67: computes `err_vec` using `err_vec=A*x-y`.
- Lines 70: computes `g` using `g=2*(A_ct*err_vec)`.
- Lines 73: computes `x_prox` using `x_prox=soft_thr(x-g/L,thr)`.
- Lines 76: computes `t_new` using `t_new=0.5*(1+sqrt(1+4*t^2))`.
- Lines 82: computes `soln_norm` using `soln_norm=norm(x_prox,2)`.
- Lines 83: computes `step_norm` using `step_norm=norm(x_prox-x_old,2)`.
- Lines 84: computes `converged` using `converged=(step_norm<step_norm_tol*soln_norm)`.
- Lines 88: computes `iter_count` using `iter_count=iter_count+1`.
- Lines 105: computes `thr_upper` using `thr_upper=thr`.

### Local helper functions

- Line 144: `grumble()` — `function grumble(A,y,nnzt)`.
  - Representative operation: `if (~isnumeric(A))||(~ismatrix(A))`.
  - Representative operation: `error('A must be a numeric matrix.')`.

## Parameters / inputs

- A -a real or complex matrix
- y -a real or complex column vector
- nnzt -the target for the number of
- non-zeroes in the solution

## Outputs

- x -a real or complex vector
- err -squared 2-norm of the fitting
- error divided by the squared
- 2-norm of the solution
- reg -1-norm of the solution

## Implementation structure

- L1 norm Tikhonov regularised solver for A*x=y where
- A is an ill-conditioned matrix. The error functional
- is norm(A*x-y,2)^2+lambda*norm(x,1), it is minimised
- using the FISTA algorithm. The user specifies the de-
- sired number of non-zeroes, lambda parameter is then
- found by bracketing / bisection. Syntax:
- [x,err,reg]=tikhol1n(A,y,nnzt)
- A -a real or complex matrix
- y -a real or complex column vector
- nnzt -the target for the number of
- non-zeroes in the solution
- x -a real or complex vector

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ctranspose()`, `normest()`, `sign()`, `int2str()`, `soft_thr()`, `nnz()`, `ismember()`, `num2str()`, `ismatrix()`, `iscolumn()`, `isscalar()`.
