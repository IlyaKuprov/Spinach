# examples/fundamentals/quadratures/product_quadratures_2.m

- MATLAB implementation: [examples/fundamentals/quadratures/product_quadratures_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/quadratures/product_quadratures_2.m)

## Purpose

Compare product-quadrature, Lie-group, and Runge-Kutta-Munthe-Kaas (RKMK) integrators for a chirped oscillator with radiation damping. The generator depends on both time and the current magnetisation.

## Model and numerical question

The Bloch-Maxwell state is a three-component magnetisation vector. The source sets the chirp-rate parameter to `2*pi*400` (commented as 400 Hz/s), longitudinal and transverse relaxation parameters to 10 Hz each, and radiation damping to 40 Hz. The initial vector is the z-axis magnetisation rotated by 178 degrees using `euler2dcm(0, pi*178/180, 0)`. The integration interval is 0.5 s.

The question is how each method's final magnetisation differs from a high-resolution RKMK-DP8 propagation as the number of time-grid points changes. The radiation-damping term is attributed by the source to Bloembergen and Pound ([DOI: 10.1103/PhysRev.95.8](https://doi.org/10.1103/PhysRev.95.8)).

## Method and checks

The reference stores the trajectory at 4,096 equally spaced points and uses RKMK-DP8. Ten benchmark grid sizes are formed by rounding `2^x` for ten equally spaced values of `x` from 8 to 11. The methods are piecewise-constant left-edge (`PWCL`), Lie-group orders 2 and 4 (`LG2` and `LG4`), the alternate fourth-order `LG4A` formula attributed in the source to Casas-Iserles Appendix A1, and `RKMK4`, `RKMK-DP5`, and `RKMK-DP8`.

For each grid and method the endpoint error is `norm(mu - mu_ref) / norm(mu_ref)`. The source plots the reference magnetisation components and the error curves on log-log axes, with the displayed error range limited to `1e-7` through `1`. It also fits a straight line to `log(error)` versus `log(grid size)` using benchmark entries 3 through 7 and prints the negative fitted slope as an empirical order. The plot bounds are not acceptance thresholds, and the fitted value is a diagnostic rather than a pass/fail test.

## Output and limitations

The example produces a trajectory plot, a relative-error comparison, and fitted empirical slopes. No measured error values or slopes are asserted here. The reference is a numerical calculation with RKMK-DP8, which is also one of the benchmarked methods; it is not an independent analytic solution. The source cites the radiation-damping model but does not provide a separate quantitative validation criterion.

[Source example](../../../../../../examples/fundamentals/quadratures/product_quadratures_2.m).
