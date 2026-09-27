# examples/benchmarks/iserstep_bench_hiord.m

## Status

This page documents a historical Spinach example, not a current runnable example. The source file `examples/benchmarks/iserstep_bench_hiord.m` was deleted in commit `c70f9b30` (“Deleting a duplicate example file”) and is absent from the current checkout. The historical function signature was `iserstep_bench_hiord()`; it is not an available entry point in the current source tree.

## Historical purpose and physics

The example compared higher-order `iserstep` methods for a chirped-frequency oscillator with radiation damping and a generator depending on both time and magnetisation. The comments cite Bloembergen and Pound for radiation damping: [Phys. Rev. 95, 8 (1954)](https://doi.org/10.1103/PhysRev.95.8).

In the historical source, the chirp rate was `2*pi*400`, the longitudinal and transverse relaxation rates were both `10`, and the radiation-damping rate was `40`. Its Bloch–Maxwell generator combined a chirped transverse precession/relaxation matrix with a magnetisation-dependent radiation-damping matrix, with the `-1i` Liouvillian factors included. The initial magnetisation was rotated from the positive z direction by 178 degrees. These values and the setup describe the old example only; they are not a current benchmark configuration.

## Historical method and outputs

A 4096-point RKMK-DP8 propagation over 0.5 seconds supplied the reference final magnetisation. The source then compared terminal-state relative errors on ten grids from `ceil(2.^linspace(8,10.5,10))`: PWCL, LG2, LG4, and LG4A through `iserstep`, and RKMK4, RKMK-DP5, and RKMK-DP8 through `step`. It plotted the reference magnetisation trajectory and relative error versus grid size, and printed empirical convergence orders by fitting log(error) against log(grid size) over the last third of the grids.

Those are outputs and comparisons described by the historical source, not instructions or confirmation that the deleted example can be run in the current checkout. The historical source is available in repository history as `c70f9b30^:examples/benchmarks/iserstep_bench_hiord.m`.
