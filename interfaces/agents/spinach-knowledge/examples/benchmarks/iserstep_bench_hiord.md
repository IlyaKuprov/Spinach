# examples/benchmarks/iserstep_bench_hiord.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/benchmarks/iserstep_bench_hiord.m`
- Signature: `iserstep_bench_hiord()`
- Total lines: 191

## Purpose

Benchmarks iserstep higher-order methods on a chirped-frequency oscillator with radiation damping, that has a state-dependent, and a time-dependent evolution generator. Syntax: iserstep_bench_hiord()

## Physical / mathematical content

- Benchmark examples. These files stress-test Spinach performance, scaling, and numerical throughput on representative spin-dynamics workloads, so runtime, memory pressure, and solver/pathway choices are part of the intended content.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file also defines local helper function(s): `run_method()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- (none) -produces a set of diagnostic plots, and prints the empirical
- convergence orders for each method

## Implementation structure

- Benchmarks iserstep higher-order methods on a chirped-frequency oscillator
- with radiation damping, that has a state-dependent, and a time-dependent
- evolution generator. Syntax:
- iserstep_bench_hiord()
- (none) -produces a set of diagnostic plots, and prints the empirical
- convergence orders for each method
- Set the chirp rate
- Set the relaxation rates
- Set the radiation damping rate
- Bootstrap the object
- Make Bloch-Maxwell generator (Liouvillian, including -1i factors)
- Set the initial magnetisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `bootstrap()`, `euler2dcm()`, `mu_traj()`, `step()`, `run_method()`, `bench()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `kylabel()`, `klegend()`, `ylim()`, `set()`, `orders()`, `polyfit()`.
