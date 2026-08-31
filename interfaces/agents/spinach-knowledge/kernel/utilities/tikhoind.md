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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 36-37: Check consistency; implemented by `grumble(K,D,y,lam)`.
- Lines 39-40: Analytical solution; implemented by `x=((K'*K)+lam*(D'*D))\(K'*y)`.
- Lines 42-43: Error and regularisation signals; implemented by `if nargout>1, err=norm(K*x-y,2)^2; end`.

### Control flow inferred from the code

- Line 43: conditional branch on `nargout>1, err=norm(K*x-y,2)^2; end`.
- Line 44: conditional branch on `nargout>2, reg=norm(D*x,2)^2; end`.

### Key state/data transformations

- Lines 40: computes `x` using `x=((K'*K)+lam*(D'*D))\(K'*y)`.

### Local helper functions

- Line 49: `grumble()` — `function grumble(K,D,y,lam)`.
  - Representative operation: `if (~isnumeric(K))||(~isnumeric(D))|| (~isnumeric(y))||(~isnumeric(lam))`.
  - Representative operation: `(~isnumeric(y))||(~isnumeric(lam))`.

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
