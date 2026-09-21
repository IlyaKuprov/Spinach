# kernel/utilities/tikhoind.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/tikhoind.m`
- Signature: `[x,err,reg]=tikhoind(K,D,y,lam)`
- Total lines: 66

## Purpose

Analytical Tikhonov regularised solution to K*x=y without any constraints (sign-indefinite output). Syntax: [x,err,reg]=tikhoind(K,D,y,lam)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- K -kernel matrix, may be complex, may be non-square
- D -regularisation matrix
- y -a column vector, may be complex
- lam -Tikhonov regularisation parameter

## Outputs

- x -a real vector, a minimum of
- norm(K*x-y,2)^2+lambda*norm(D*x,2)^2
- err -error signal norm(K*x-y,2)^2
- reg -regularisation signal norm(D*x,2)^2
- Note: for best numerical performance, scale K to have approxima-
- tely unit 2-norm, and y to have approximately unit 1-norm.
- Note: see tikhonov.m for the positive-constraned solver.

## Implementation structure

- Analytical Tikhonov regularised solution to K*x=y without any
- constraints (sign-indefinite output). Syntax:
- [x,err,reg]=tikhoind(K,D,y,lam)
- K -kernel matrix, may be complex, may be non-square
- D -regularisation matrix
- y -a column vector, may be complex
- lam -Tikhonov regularisation parameter
- x -a real vector, a minimum of
- norm(K*x-y,2)^2+lambda*norm(D*x,2)^2
- err -error signal norm(K*x-y,2)^2
- reg -regularisation signal norm(D*x,2)^2
- Note: for best numerical performance, scale K to have approxima-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`.
