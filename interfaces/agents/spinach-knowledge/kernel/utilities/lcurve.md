# kernel/utilities/lcurve.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/lcurve.m`
- Signature: `lam_opt=lcurve(lam,err,reg,mode)`
- Total lines: 165

## Purpose

L-curve analysis function. Syntax: lam_opt=lcurve(lam,err,reg,mode)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- lam -row vector of regularisation parameters, must
- be positive and in ascending order
- err -row vector of least squares errors, must be
- positive and increasing with lam
- reg -row vector of regularisation functional values,
- must be positive and decreasing with lam. This
- is the regularisation functional itself, not the
- penalty term of the error functional: when the
- optimiser reports lam*||L*x||^2, divide it by lam
- once before calling this function
- mode -'log' for logarithmic coordinates and 'linear'
- for linear ones; 'log' is recommended

## Outputs

- lam_opt -the regularisation parameter at the point
- of the maximum curvature of the L-curve
- Notes: the corner is the point of the greatest curvature, which for a
- smooth asymmetric bend is not the intersection of the asymptotes;
- the criterion locates the regularisation parameter to within a
- factor of a few and is not convergent in the zero noise limit
- (Vogel, SIAM J. Numer. Anal. 34, 1996). Treat the answer as an
- order of magnitude estimate and inspect the plotted curve.
- The curvature maximum must fall inside the sampled interval; if
- it falls on either end, the corner is outside the range and this
- function refuses to return the endpoint as an answer.
- This function requires the Curve Fitting Toolbox.

## Implementation structure

- L-curve analysis function. Syntax:
- lam_opt=lcurve(lam,err,reg,mode)
- lam -row vector of regularisation parameters, must
- be positive and in ascending order
- err -row vector of least squares errors, must be
- positive and increasing with lam
- reg -row vector of regularisation functional values,
- must be positive and decreasing with lam. This
- is the regularisation functional itself, not the
- penalty term of the error functional: when the
- optimiser reports lam*||L*x||^2, divide it by lam
- once before calling this function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `any()`, `diff()`, `log10()`, `spapi()`, `optknt()`, `fnval()`, `subplot()`, `kxlabel()`, `kylabel()`, `fdvec()`, `set()`, `kappa()`, `lam()`, `err()`, `reg()`.
