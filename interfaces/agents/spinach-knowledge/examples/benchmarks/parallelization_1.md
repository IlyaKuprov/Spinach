# examples/benchmarks/parallelization_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/benchmarks/parallelization_1.m`
- Signature: `parallelization_1()`
- Total lines: 67

## Purpose

Parallelization test: multi-threaded evaluation of observables in Hilbert space time propagation. For further information, see: http://dx.doi.org/10.1063/1.3679656 Spin system of 3-phenylmethylene-1H,3H-naphtho-[1,8-c,d]-pyran-1-one. Source: Penchav, et al., Spec. Acta Part A, 78 (2011) 559-565.

## Physical / mathematical content

- Benchmark examples. These files stress-test Spinach performance, scaling, and numerical throughput on representative spin-dynamics workloads, so runtime, memory pressure, and solver/pathway choices are part of the intended content.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `ncores()`, `feature()`, `delete()`, `gcp()`, `parpool()`, `pause()`, `evolution()`, `num2str()`.
