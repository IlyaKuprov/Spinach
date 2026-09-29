# examples/fundamentals/quadratures/mas_benchmark.m

- MATLAB implementation: [examples/fundamentals/quadratures/mas_benchmark.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/quadratures/mas_benchmark.m)

## Purpose

Compare four time-stepping schemes for one period of a magic-angle-spinning rotor, using a finely sampled propagation as the reference. The source describes the calculation as a benchmark; it does not define a pass/fail test.

## System and numerical question

The model is a two-proton system at 14.1 T, with scalar Zeeman values 5.0 and -2.0 and coordinates `0 0 0` and `0 3.9 0.1` as configured in the source. It uses the `sphten-liouv` basis, no basis approximation, projection `{+1}`, and the NMR assumption. The initial state is proton `L+`. The rotor speed is 50,000 Hz, so one period is 20 microseconds; the rotor orientation is sampled at the magic angle `atan(sqrt(2))`.

The question is how the relative final-state error changes with rotor-grid size for left-point and midpoint piecewise-constant propagation, and second- and fourth-order Lie-group propagation.

## Method and checks

The reference uses 8,193 equally spaced rotor points (8,192 intervals). It advances through the stored Hamiltonian samples in two-interval steps with a three-sample `step` call. Benchmark grids contain `2^p + 1` points for `p = 4,...,12`. Each method propagates the same initial state over one period. For each grid, the plotted quantity is `norm(rho_ref - rho_method) / norm(rho_ref)`.

The source plots error against grid-point count on logarithmic axes and clips the displayed vertical range to `1e-13` through `1e-2`. These are plot limits, not acceptance thresholds; the example contains no pass/fail criterion.

## Output and limitations

The output is a log-log comparison plot for the four schemes. No numerical error values or measured convergence rates are stated by the source, and this page does not assert them. The reference is a fine-grid propagation using the same `step` machinery, not an analytic solution. The source header estimates runtime as seconds; this is not a measured runtime here.

[Source example](../../../../../../examples/fundamentals/quadratures/mas_benchmark.m).
