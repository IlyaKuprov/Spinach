# examples/benchmarks/parallelization_2.m

- Signature: `parallelization_2()`

## Purpose

Parallelization test: multi-threaded evaluation of observables in Hilbert-space time propagation for a pyrene radical spin system at low field. For further information, see: http://dx.doi.org/10.1063/1.3679656

## Physical / mathematical content

- Models pyrene radical at 50 µT with one proton and two electron spins.

## Numerical / algorithmic content

- Uses the Zeeman Hilbert-space formalism and lab-frame Hamiltonian; times a 200-step observable propagation while varying the parallel-pool size up to the available core count.

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
