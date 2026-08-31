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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Read the spin system properties (vacuum DFT calculation); implemented by `options.no_xyz=1`.
- Lines 16-17: Magnet field; implemented by `sys.magnet=50e-6`.
- Lines 19-20: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 23-24: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 27-28: Assumptions; implemented by `spin_system=assume(spin_system,'labframe')`.
- Lines 30-31: Hamiltonian operator; implemented by `[H,Q]=hamiltonian(spin_system)`.
- Lines 34-35: Initial state; implemented by `rho=operator(spin_system,'Lz','E')`.
- Lines 37-38: Parallel propagation benchmark, 200 steps; implemented by `ncores=[2 4 8 16 32 64 128 256 512 1024]`.

### Control flow inferred from the code

- Line 40: `for` loop over `n=ncores`.

### Key state/data transformations

- Lines 13: computes `options.no_xyz` using `options.no_xyz=1`.
- Lines 14-15: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/pyrene_cation.log'), {{'E','E'},{'H','1H'}},[0 0],options)`.
- Lines 17: computes `sys.magnet` using `sys.magnet=50e-6`.
- Lines 20: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `[H,Q]` using `[H,Q]=hamiltonian(spin_system)`.
- Lines 32: computes `H` using `H=H+orientation(Q,[pi/3,pi/4,pi/5])`.
- Lines 35: computes `rho` using `rho=operator(spin_system,'Lz','E')`.
- Lines 38: computes `ncores` using `ncores=[2 4 8 16 32 64 128 256 512 1024]`.

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
