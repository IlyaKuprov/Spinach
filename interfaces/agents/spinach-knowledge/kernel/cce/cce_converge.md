# kernel/cce/cce_converge.m

## Native applicability

This retained expert convergence utility has no native `clexpn` invocation contract in the current implementation. Do not pass native experiment parameters to it or infer a registered experiment from its two-argument signature.

For native simulations, call `clexpn(spin_system,pulse_sequence,parameters,assumptions)` explicitly for each refinement and compare numerical answers together with `diagnostics.valid` and `diagnostics.convergence`. The context's omitted-order policy is not an independent bath-radius, graph-threshold or sample-count convergence certificate.

## Convergence meaning

Refine cluster order, bath extent, graph cutoff, bath-state samples and time resolution independently while holding other physical choices fixed. Compare time grids only on declared common samples. Use explicit absolute/relative observable tolerances and a bounded evaluation budget. Several successive levels agreeing establishes the requested ladder criterion, not universal physical convergence. Bath-configuration convergence requires independently prepared physical realizations rather than a fabricated configuration generator.

## Source

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_converge.m)
