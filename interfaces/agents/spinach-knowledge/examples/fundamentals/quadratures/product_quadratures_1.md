# examples/fundamentals/quadratures/product_quadratures_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/quadratures/product_quadratures_1.m`
- Signature: `product_quadratures_1()`
- Total lines: 135

## Purpose

Accuracy test for Lie-group product quadratures as a function of discretisation step in the E1000B Veshtort-Griffin pulse. Calculation time: seconds

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnetic field; implemented by `sys.magnet=14.1`.
- Lines 15-18: Isotopes; implemented by `sys.isotopes={ '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1…`.
- Lines 20-21: Zeeman interactions; implemented by `inter.zeeman.scalar=num2cell(linspace(-4,4,31))`.
- Lines 23-24: Couplings; implemented by `inter.coupling.scalar=cell(31)`.
- Lines 29-30: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 42-43: Hamiltonian superoperator; implemented by `H=hamiltonian(spin_system)`.
- Lines 45-46: Control operator; implemented by `Lx=operator(spin_system,'Lx','1H')`.
- Lines 48-49: Pulse duration; implemented by `duration=1e-2`.
- Lines 51-52: Accurate reference calculation; implemented by `rho_ref=state(spin_system,'Lz','1H')`.
- Lines 61-62: Plotting; implemented by `kfigure(); scale_figure([1.50 0.60])`.
- Lines 69-70: Benchmark arrays; implemented by `npi=100:100:1000; bench=zeros(numel(npi),4)`.
- Lines 72-73: Benchmarking loop; implemented by `parfor k=1:numel(npi)`.
- Lines 75-76: Left point integrator; implemented by `rho_left=state(spin_system,'Lz','1H')`.
- Lines 83-84: Midpoint integrator; implemented by `rho_mid=state(spin_system,'Lz','1H')`.
- Lines 92-93: Two-point integrator; implemented by `rho_two=state(spin_system,'Lz','1H')`.
- Lines 101-102: Three-point integrator; implemented by `rho_thr=state(spin_system,'Lz','1H')`.

### Control flow inferred from the code

- Line 25: `for` loop over `n=1:30`.
- Line 55: `for` loop over `n=1:2:(numel(amps)-2)`.
- Line 73: `parfor` loop over `k=1:numel(npi)`.
- Line 79: `for` loop over `n=1:(np-1)`.
- Line 87: `for` loop over `n=1:(np-1)`.
- Line 96: `for` loop over `n=1:(np-1)`.
- Line 105: `for` loop over `n=1:2:(numel(amps)-2)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 16-18: computes `sys.isotopes` using `sys.isotopes={ '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1…`.
- Lines 21: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=num2cell(linspace(-4,4,31))`.
- Lines 24: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(31)`.
- Lines 26: computes `inter.coupling.scalar{n,n+1}` using `inter.coupling.scalar{n,n+1}=10`.
- Lines 30: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 31: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 32: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 33: computes `bas.space_level` using `bas.space_level=1`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 43: computes `H` using `H=hamiltonian(spin_system)`.
- Lines 46: computes `Lx` using `Lx=operator(spin_system,'Lx','1H')`.
- Lines 49: computes `duration` using `duration=1e-2`.
- Lines 52: computes `rho_ref` using `rho_ref=state(spin_system,'Lz','1H')`.
- Lines 53: computes `np` using `np=2000; dt=duration/(np-1)`.
- Lines 54: computes `amps` using `amps=vg_pulse('E1000B',2*np-1,duration)`.
- Lines 63: computes `time_axis` using `time_axis=linspace(0,duration,2*np-1)'`.
- Lines 70: computes `npi` using `npi=100:100:1000; bench=zeros(numel(npi),4)`.

## Implementation structure

- Accuracy test for Lie-group product quadratures as a function of
- discretisation step in the E1000B Veshtort-Griffin pulse.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping
- Assumptions
- Hamiltonian superoperator
- Control operator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `state()`, `vg_pulse()`, `step()`, `amps()`, `kfigure()`, `scale_figure()`, `subplot()`, `kylabel()`, `kxlabel()`, `npi()`.
