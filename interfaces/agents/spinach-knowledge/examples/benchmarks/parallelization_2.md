# examples/benchmarks/parallelization_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/benchmarks/parallelization_2.m`
- Signature: `parallelization_2()`
- Total lines: 47

## Purpose

Parallelization test: multi-threaded evaluation of observables in Hilbert space time propagation for pyrene radical spin system at low field. For further information, see:

## Physical / mathematical content

- Benchmark examples. These files stress-test Spinach performance, scaling, and numerical throughput on representative spin-dynamics workloads, so runtime, memory pressure, and solver/pathway choices are part of the intended content.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Parallelization test: multi-threaded evaluation of observables in Hilbert
- space time propagation for pyrene radical spin system at low field. For
- further information, see:
- Read the spin system properties (vacuum DFT calculation)
- Magnet field
- Basis set
- Spinach housekeeping
- Assumptions
- Hamiltonian operator
- Initial state
- Parallel propagation benchmark, 200 steps

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `orientation()`, `operator()`, `ncores()`, `feature()`, `delete()`, `gcp()`, `parpool()`, `pause()`, `evolution()`, `num2str()`.
