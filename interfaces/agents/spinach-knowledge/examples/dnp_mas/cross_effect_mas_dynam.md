# examples/dnp_mas/cross_effect_mas_dynam.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_mas/cross_effect_mas_dynam.m`
- Signature: `cross_effect_mas_dynam()`
- Total lines: 94

## Purpose

A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different): Spin system trajectory analysis within the first rotor period for a single crystal. Calculation time: seconds

## Physical / mathematical content

- MAS DNP examples. These files model microwave-driven electron-nuclear polarisation transfer under magic-angle spinning, combining rotor-synchronised anisotropic interactions, relaxation, microwave irradiation, and powder/rotor averaging.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Magnet field; implemented by `sys.magnet=9.394`.
- Lines 18-19: Spin specification; implemented by `sys.isotopes={'E','E','1H'}`.
- Lines 21-22: Interactions; implemented by `inter.zeeman.eigs{1}=[2.0094 2.0060 2.0017]`.
- Lines 35-36: Relaxation parameters; implemented by `inter.relaxation={'nottingham'}`.
- Lines 45-46: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 49-50: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 53-54: Stack generation parameters; implemented by `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 62-63: Rotor stack generation; implemented by `H=rotor_stack(spin_system,parameters,'esr')`.
- Lines 65-66: Microwave operator; implemented by `Hmw=operator(spin_system,'Lx','E')`.
- Lines 68-69: Zeeman offset operator; implemented by `HzE=operator(spin_system,'Lz','E')`.
- Lines 71-72: Experiment parameters; implemented by `parameters.mw_pwr=0.85e6`.
- Lines 76-77: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 79-80: Thermal equilibrium state; implemented by `rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 82-83: Rotor period trajectory; implemented by `nsteps=numel(H); rho=zeros(numel(rho_eq),nsteps,'like',1i); rho(:,1)=rho_eq`.
- Lines 90-91: Trajectory analysis; implemented by `kfigure(); trajan(spin_system,rho,'level_populations')`.

### Control flow inferred from the code

- Line 84: `for` loop over `n=2:nsteps`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=9.394`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'E','E','1H'}`.
- Lines 22: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.0094 2.0060 2.0017]`.
- Lines 23: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[0.00 0.00 0.00]`.
- Lines 24: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[2.0094 2.0060 2.0017]`.
- Lines 25: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=pi*[107 108 124]/180`.
- Lines 26: computes `inter.zeeman.eigs{3}` using `inter.zeeman.eigs{3}=[0.00 0.00 0.00]`.
- Lines 27: computes `inter.zeeman.euler{3}` using `inter.zeeman.euler{3}=[0.00 0.00 0.00]`.
- Lines 28: computes `inter.coupling.eigs` using `inter.coupling.eigs=cell(3,3)`.
- Lines 29: computes `inter.coupling.euler` using `inter.coupling.euler=cell(3,3)`.
- Lines 30: computes `inter.coupling.eigs{1,2}` using `inter.coupling.eigs{1,2}=[23.0e6 -11.5e6 -11.5e6]`.
- Lines 31: computes `inter.coupling.euler{1,2}` using `inter.coupling.euler{1,2}=pi*[0.00 135 0.00]/180`.
- Lines 32: computes `inter.coupling.eigs{1,3}` using `inter.coupling.eigs{1,3}=[1.5e6 -0.75e6 -0.75e6]`.
- Lines 33: computes `inter.coupling.euler{1,3}` using `inter.coupling.euler{1,3}=[0.00 0.00 0.00]`.
- Lines 36: computes `inter.relaxation` using `inter.relaxation={'nottingham'}`.
- Lines 37: computes `inter.nott_r1e` using `inter.nott_r1e=1/0.3e-3`.
- Lines 38: computes `inter.nott_r1n` using `inter.nott_r1n=1/4.0`.
- Lines 39: computes `inter.nott_r2e` using `inter.nott_r2e=1/1.0e-6`.

## Implementation structure

- A MAS DNP simulation performed as described in Fred Mentink-
- Vigier's paper (Spinach rotation conventions are different):
- Spin system trajectory analysis within the first rotor period
- for a single crystal.
- Calculation time: seconds
- Magnet field
- Spin specification
- Interactions
- Relaxation parameters
- Basis set
- Spinach housekeeping
- Stack generation parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rotor_stack()`, `operator()`, `relaxation()`, `equilibrium()`, `hamiltonian()`, `assume()`, `rho()`, `step()`, `kfigure()`, `trajan()`.
