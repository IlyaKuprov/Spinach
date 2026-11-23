# kernel/cce/cce_converge.m

## Native convergence contract

Use `[answer,diagnostics]=cce_converge(spin_system,parameters,pulse_sequence,assumptions)`. The experiment and assumptions are explicit, exactly as in clexpn; the primary result is the unchanged numeric baseline, never a wrapper or a silently substituted refined result.

Declare independent ladders in parameters.cce.convergence using the paths cce.order, cce.bath_spins, bath_samples, or times. Every ladder starts at the supplied baseline. Every refinement is compared with that baseline, including its numerical masks and error bounds. Orders above the finite bath size evaluate the same full-bath family. Time ladders retain every preceding time exactly; interpolation is not used.

Supply atol, rtol, and max_work. Sampling ladders require a held-out bath-state seed in checks; nested cce overrides replace only the specified algorithm controls. Shared checkpoints are rejected before any evaluation. The second output contains convergence (converged, limit_reached, or not_converged), the original input_parameters, and convergence_report with comparisons, work, coverage, and budgets. The diagnostic reduction and ensemble remain those of the unchanged baseline.

## Convergence meaning

Refine cluster order, bath extent, graph cutoff, bath-state samples, and time resolution independently while holding other physical choices fixed. Compare time grids only on declared common samples. Use explicit absolute/relative observable tolerances and a bounded evaluation budget. Several successive levels agreeing establishes the requested ladder criterion, not universal physical convergence. Bath-configuration convergence requires independently prepared physical realisations rather than a fabricated configuration generator.

## Source

[MATLAB source](../../../../../kernel/cce/cce_converge.m)

Order-ladder entries must be numeric, real, finite nonnegative integer scalars. They are validated before finite-bath saturation, which must not turn NaN, infinity, fractional values, or text into a valid order.

Order saturation uses the resolved baseline bath IDs, including the default complement of observed spins. It does not require an explicit cce.bath_spins field or insert one into the recorded caller input.

The maximum sampling uncertainty across shared times includes missing estimates: if any requested-time contribution is `NaN`, the reported `sampling_error` remains unavailable. A finite estimate at another time cannot hide it. Coverage and numerical-validity checks remain independent and unchanged.

The experiment function handle and explicit assumption arguments are checked locally before any baseline or refinement is evaluated.
