# kernel/utilities/lcurve.m

- Signature: `lam_opt=lcurve(lam,err,reg,mode)`

## Purpose

L-curve analysis function. Syntax: lam_opt=lcurve(lam,err,reg,mode)

## Physical / mathematical content

- Plots the trade-off between least-squares error and regularisation functional and estimates the parameter at the L-curve maximum-curvature corner.

## Numerical / algorithmic content

- Resamples quintic splines of the log10 error and regularisation values at 1000 points, then estimates curvature using finite differences in the requested `log` or `linear` coordinates. The maximum is searched away from stencil edges and must lie inside the sampled range.

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

- Validates equal-length positive row vectors, at least six increasing `lam` values, and `mode` (`log` or `linear`); non-monotone error or regularisation values trigger a warning.
- Plots the L-curve and curvature, returns the parameter at the interior curvature maximum, and errors when the corner lies near a sampled-range edge.
- Requires MATLAB Curve Fitting Toolbox.
