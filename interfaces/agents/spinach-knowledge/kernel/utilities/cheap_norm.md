# kernel/utilities/cheap_norm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/cheap_norm.m`
- Signature: `n=cheap_norm(A,t,itmax)`
- Total lines: 182

## Purpose

The cheapest norm for various representations of matrices. CUDA stores matrices by rows, Matlab by columns, and polyadic objects can only multiply vectors, in which case Algorithm 2.4 from Hig- ham and Tisseur's paper: is used. Syntax: n=cheap_norm(A,t,itmax)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- A -a matrix, or a polyadic representation thereof
- t -(optional) number of probe columns in the poly-
- adic norm estimator, defaults to 1
- itmax -(optional) maximum number of estimator iterati-
- ons, defaults to 5

## Outputs

- n -infinity-norm for GPU arrays, 1-norm for CPU arrays,
- and a lower-bound 1-norm estimate for polyadics
- Note: some norms are vastly more expensive than others, this
- function uses the cheapest ones available.

## Implementation structure

- The cheapest norm for various representations of matrices. CUDA
- stores matrices by rows, Matlab by columns, and polyadic objects
- can only multiply vectors, in which case Algorithm 2.4 from Hig-
- ham and Tisseur's paper:
- is used. Syntax:
- n=cheap_norm(A,t,itmax)
- A -a matrix, or a polyadic representation thereof
- t -(optional) number of probe columns in the poly-
- adic norm estimator, defaults to 1
- itmax -(optional) maximum number of estimator iterati-
- ons, defaults to 5
- n -infinity-norm for GPU arrays, 1-norm for CPU arrays,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `randi()`, `any()`, `idx()`, `sign()`, `all()`, `row_scores()`, `ismember()`, `isscalar()`.
