# kernel/utilities/lcurve.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/lcurve.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/lcurve.m)

## Purpose

L-curve analysis function that locates the optimal regularisation parameter at the point of maximum curvature of the L-curve, given a sweep of regularisation parameters with corresponding least squares errors and regularisation functional values.

## Behaviour

- Syntax: `lam_opt=lcurve(lam,err,reg,mode)`.
- Input consistency is enforced by an internal `grumble` subfunction: `lam`, `err` and `reg` must be row vectors of positive, real, finite numbers of equal length, with at least six elements; `lam` must be in ascending order; `mode` must be `'log'` or `'linear'`.
- If `err` does not increase or `reg` does not decrease with `lam`, a warning is issued ("err should increase and reg should decrease with lam, inspect the sweep.") but execution continues.
- All three vectors are converted to base-10 logarithms and resampled onto 1000 points over the sampled interval using quartic (order-five) splines (`spapi` with `optknt(...,5)` knots), then converted back to linear coordinates.
- A two-panel figure is produced: the L-curve (`err` versus `reg`) on the left and the curvature versus `lam` on the right, both with log-scaled axes; the located optimum is marked with a red circle on both panels.
- Derivatives are computed with `fdvec` using 5-point stencils, either in logarithmic coordinates (`'log'` mode) or linear coordinates (`'linear'` mode); any other mode string raises an error.
- The signed curvature is computed as `kappa=(xp.*ypp-yp.*xpp)./((xp.^2+yp.^2).^(3/2))`.
- The maximum curvature is searched only away from the stencil ends, with a margin of `ceil(numel(kappa)/40)` points excluded at each end.
- If the curvature maximum falls within `2*margin` points of either end, the function errors out with "L-curve corner is outside the regularisation parameter range, widen it." rather than returning an endpoint.
- Notes from the header: the corner is the point of greatest curvature, which for a smooth asymmetric bend is not the intersection of the asymptotes; the criterion locates the regularisation parameter to within a factor of a few and is not convergent in the zero noise limit; the answer should be treated as an order of magnitude estimate and the plotted curve inspected.
- The function requires the Curve Fitting Toolbox.

## Inputs and outputs

**Inputs:**

- `lam` — row vector of regularisation parameters, positive and in ascending order.
- `err` — row vector of least squares errors, positive and increasing with `lam`.
- `reg` — row vector of regularisation functional values, positive and decreasing with `lam`. This is the regularisation functional itself, not the penalty term of the error functional: when the optimiser reports `lam*||L*x||^2`, it must be divided by `lam` once before calling this function.
- `mode` — `'log'` for logarithmic coordinates and `'linear'` for linear ones; `'log'` is recommended.

**Outputs:**

- `lam_opt` — the regularisation parameter at the point of maximum curvature of the L-curve.

## References

- Vogel, C. R., *SIAM Journal on Numerical Analysis*, **34**, 1996 (cited in the source header regarding non-convergence of the L-curve criterion in the zero noise limit).
- Spinach Wiki: [lcurve.m](https://spindynamics.org/wiki/index.php?title=lcurve.m)
