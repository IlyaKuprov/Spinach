# examples/fundamentals/quadratures/grid_leb_vs_etc.m

## Purpose and question

This example compares spherical-harmonic integration-error profiles for a Lebedev grid and six other grid constructions. The source comment calls it a heuristic-versus-Lebedev bake-off, but a source comment is not evidence that one family wins; the plotted data must be produced to assess that question.

## Benchmark configuration

Every profile is requested from `grid_test` for spherical-harmonic ranks `4:2:60` using `Y_lm`. The compared configurations are:

- Lebedev rank 29, loaded with its angles and weights; the legend labels it 302 points.
- Repulsion with the same number of points as that Lebedev grid, generated with parameter 3 and 10000 iterations. Voronoi weights are computed from the unit-sphere coordinates derived from beta and gamma, then divided by `4*pi`. Its legend also says 302 points.
- ZCWn from `grid_fibon('zcwn',302)`, labelled 302 points.
- Igloo from `grid_igloo(17)`, labelled 328 points.
- Stoll from `grid_trian('stoll',9)`, labelled 326 points.
- ASG and SOPHE from `grid_trian` at level 9, each labelled 326 points.

The non-Lebedev grid helpers provide angles and weights to the test; the source explicitly constructs and normalises Voronoi weights for the repulsion case. The plot uses logarithmic integration-error axes, limits the displayed error range to `10^-16` through 1, and displays spherical rank from 0 to 64.

## Output and interpretation

The function draws seven error profiles against spherical-harmonic rank. It contains no pass/fail threshold and does not print tabulated error values, so the source alone does not establish a performance ranking or quantitative error bound.

There is a label-to-construction clarification: the second plotted series is computed from the Stoll grid call, but its legend text is `EasySpin, 326 pts`. The source does not show an EasySpin call for that series. The source also does not include citations or a saved benchmark result. The plotted profiles depend on the local grid files, helper implementations, and a MATLAB run.

Source: [grid_leb_vs_etc.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/quadratures/grid_leb_vs_etc.m).
