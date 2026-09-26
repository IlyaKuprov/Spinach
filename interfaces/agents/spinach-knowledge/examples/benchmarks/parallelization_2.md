# examples/benchmarks/parallelization_2.m

- Signature: `parallelization_2()`

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
