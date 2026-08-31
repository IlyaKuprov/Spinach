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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnetic induction; implemented by `sys.magnet=14.095`.
- Lines 17-19: Spin system; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H'}`.
- Lines 21-23: Chemical shifts; implemented by `inter.zeeman.scalar={8.345,7.741,8.097,8.354,7.784,8.330, 7.059,7.941,7.466,7.326,7.466,7.941}`.
- Lines 25-26: Scalar couplings; implemented by `inter.coupling.scalar=cell(12,12)`.
- Lines 40-41: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 44-45: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 48-49: Assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 51-52: Hamiltonian operator; implemented by `H=hamiltonian(spin_system)`.
- Lines 54-55: Initial state; implemented by `rho=operator(spin_system,'Lx','all')`.
- Lines 57-58: Parallel propagation benchmark, 1000 steps; implemented by `ncores=[2 4 8 16 32 64 128 256 512 1024]`.

### Control flow inferred from the code

- Line 60: `for` loop over `n=ncores`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=14.095`.
- Lines 18-19: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H'}`.
- Lines 22-23: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={8.345,7.741,8.097,8.354,7.784,8.330, 7.059,7.941,7.466,7.326,7.466,7.941}`.
- Lines 26: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(12,12)`.
- Lines 27: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=7.8`.
- Lines 28: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=0.9`.
- Lines 29: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=7.8`.
- Lines 30: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=8.4`.
- Lines 31: computes `inter.coupling.scalar{4,6}` using `inter.coupling.scalar{4,6}=1.2`.
- Lines 32: computes `inter.coupling.scalar{5,6}` using `inter.coupling.scalar{5,6}=7.2`.
- Lines 33: computes `inter.coupling.scalar{8,9}` using `inter.coupling.scalar{8,9}=7.8`.
- Lines 34: computes `inter.coupling.scalar{8,10}` using `inter.coupling.scalar{8,10}=1.2`.
- Lines 35: computes `inter.coupling.scalar{9,10}` using `inter.coupling.scalar{9,10}=7.8`.
- Lines 36: computes `inter.coupling.scalar{10,11}` using `inter.coupling.scalar{10,11}=7.8`.
- Lines 37: computes `inter.coupling.scalar{10,12}` using `inter.coupling.scalar{10,12}=1.2`.
- Lines 38: computes `inter.coupling.scalar{11,12}` using `inter.coupling.scalar{11,12}=7.8`.
- Lines 41: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 42: computes `bas.approximation` using `bas.approximation='none'`.

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
