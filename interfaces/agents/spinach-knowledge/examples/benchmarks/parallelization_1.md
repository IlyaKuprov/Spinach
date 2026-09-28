# examples/benchmarks/parallelization_1.m

- Signature: `parallelization_1()`

## Purpose

Parallelization test: multi-threaded evaluation of observables in Hilbert space time propagation. For further information, see: http://dx.doi.org/10.1063/1.3679656 Spin system of 3-phenylmethylene-1H,3H-naphtho-[1,8-c,d]-pyran-1-one. Source: Penchav, et al., Spec. Acta Part A, 78 (2011) 559-565.

## Physical / mathematical content

- Uses a 12-proton spin system in the Zeeman Hilbert-space formalism.

## Numerical / algorithmic content

- Propagates the initial transverse magnetisation for 1,000 steps while varying the parallel-pool size up to the available core count.

## Implementation structure

- Parallelization test: multi-threaded evaluation of observables in
- Hilbert space time propagation.
- For further information, see: http://dx.doi.org/10.1063/1.3679656
- Spin system of 3-phenylmethylene-1H,3H-naphtho-[1,8-c,d]-pyran-1-one.
- Source: Penchav, et al., Spec. Acta Part A, 78 (2011) 559-565.
- Magnetic induction
- Spin system
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Assumptions
