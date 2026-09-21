# kernel/utilities/tikhonov.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/tikhonov.m`
- Signature: `[x,err,reg]=tikhonov(K,D,KtK,DtD,H,y,lambda)`
- Total lines: 119

## Purpose

Tikhonov regularised solution to K*x=y with a positivity const- raint on x using regularised Newton-Raphson method. Syntax: [x,err,reg]=tikhonov(K,D,KtK,DtD,H,y,lambda)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `lsq_err()`, `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- K -kernel matrix, may be complex, may be non-square
- D -regularisation matrix, leave empty to use finite
- difference second derivative matrix
- KtK -K'*K, for repeated calls it may be faster to pre-
- compute this quantity, leave empty otherwise
- DtD -D'*D, for repeated calls it may be faster to pre-
- compute this quantity, leave empty otherwise
- H -Tikhonov Hessian 2*real(KtK+lambda*DtD), for re-
- peated calls it may be faster to precompute this
- quantity, leave empty otherwise
- y -a column vector, may be complex
- lambda -Tikhonov regularisation parameter

## Outputs

- x -a real vector, a minimum (subject to positivity)
- of norm(K*x-y,2)^2+lambda*norm(D*x,2)^2
- err -error signal norm(K*x-y,2)^2
- reg -regularisation signal norm(D*x,2)^2
- Note: for best numerical performance, scale K to have approxima-
- tely unit 2-norm, and y to have approximately unit 1-norm.
- Note: see tikhoind.m for the indeterminate solver.

## Implementation structure

- Tikhonov regularised solution to K*x=y with a positivity const-
- raint on x using regularised Newton-Raphson method. Syntax:
- [x,err,reg]=tikhonov(K,D,KtK,DtD,H,y,lambda)
- K -kernel matrix, may be complex, may be non-square
- D -regularisation matrix, leave empty to use finite
- difference second derivative matrix
- KtK -K'*K, for repeated calls it may be faster to pre-
- compute this quantity, leave empty otherwise
- DtD -D'*D, for repeated calls it may be faster to pre-
- H -Tikhonov Hessian 2*real(KtK+lambda*DtD), for re-
- peated calls it may be faster to precompute this
- quantity, leave empty otherwise

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fdmat()`, `lsq_err()`, `optimoptions()`, `inf()`, `fmincon()`, `isscalar()`.
