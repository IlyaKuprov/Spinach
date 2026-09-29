# examples/benchmarks/iserstep_bench_hiord.m

- Signature: `iserstep_bench_hiord()`

## Problem and reference

The example compares higher-order propagators on a chirped Bloch–Maxwell oscillator with both time-dependent and state-dependent evolution. Its radiation-damping term follows Bloembergen and Pound ([Phys. Rev. 95, 8 (1954)](https://doi.org/10.1103/PhysRev.95.8)). The source sets chirp rate (2pi×400), longitudinal and transverse relaxation rates (r_1=r_2=10), and radiation-damping rate (r_{rd}=40). The generator combines the chirped transverse precession/relaxation matrix with a magnetisation-dependent radiation-damping matrix, including the source's (-i) Liouvillian factors. The initial magnetisation is a 178° rotation of the positive z direction.

A 4096-point RKMK-DP8 propagation over 0.5 s supplies a reference trajectory and terminal state. For each tested grid size, the script reports terminal-state relative error, (|mu-mu_{ref}|/|mu_{ref}|), against that reference. The reference is a high-resolution numerical comparison, not an analytic exact solution.

## Methods and diagnostics

Ten grid sizes are generated as `ceil(2.^linspace(8,10.5,10))`; the step is (0.5/(n_p-1)). The script compares PWCL, LG2, LG4 and LG4A through `iserstep`, and RKMK4, RKMK-DP5 and RKMK-DP8 through `step`. It plots the reference magnetisation trajectory and relative error against grid size. For each method it fits (log(error)) against (log(n_p)) over the last third of the grid sizes and prints the negative fitted slope as an empirical convergence order.

The source describes the setup, comparison, plots, and printed orders; it does not provide fixed results in the source comments. Do not treat the page as a performance ranking or claim that a particular method wins without running the historical code under a specified environment.
