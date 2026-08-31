# examples/dnp_mas/solid_effect_mas_steady.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_mas/solid_effect_mas_steady.m`
- Signature: `solid_effect_mas_steady()`
- Total lines: 108

## Purpose

A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different): Steady state rotor period simulation for a single crystal, computed using Newton-Raphson steady state solver. Calculation time: seconds

## Physical / mathematical content

- MAS DNP examples. These files model microwave-driven electron-nuclear polarisation transfer under magic-angle spinning, combining rotor-synchronised anisotropic interactions, relaxation, microwave irradiation, and powder/rotor averaging.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Magnet field; implemented by `sys.magnet=9.403`.
- Lines 18-19: Spin specification; implemented by `sys.isotopes={'E','1H'}`.
- Lines 21-22: Interactions; implemented by `inter.zeeman.eigs{1}=[2.00614 2.00194 2.00988]`.
- Lines 29-30: Relaxation parameters; implemented by `inter.relaxation={'weizmann'}`.
- Lines 41-42: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 45-46: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 49-50: Stack generation parameters; implemented by `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 58-59: Rotor stack generation; implemented by `H=rotor_stack(spin_system,parameters,'esr')`.
- Lines 61-62: Microwave operator; implemented by `Hmw=operator(spin_system,'Lx','E')`.
- Lines 64-65: Zeeman offset operator; implemented by `HzE=operator(spin_system,'Lz','E')`.
- Lines 67-68: Experiment parameters; implemented by `parameters.mw_pwr=0.85e6`.
- Lines 72-73: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 75-76: Rotor period integration; implemented by `nsteps=numel(H); P=speye(size(H{1}))`.
- Lines 84-85: Steady state; implemented by `rho_st=steady(spin_system,P,[],'newton')`.
- Lines 87-88: Rotor period trajectory; implemented by `nsteps=numel(H); rho=zeros(numel(rho_st),nsteps,'like',1i); rho(:,1)=rho_st`.
- Lines 95-96: Trajectory analysis; implemented by `kfigure(); trajan(spin_system,rho,'level_populations')`.
- Lines 98-99: Thermal equilibrium state; implemented by `rho_eq=equilibrium(spin_system,hamiltonian(assume(spin_system,'labframe'),'left'))`.
- Lines 101-102: Enhancement factor; implemented by `Hz_dnp=mean(state(spin_system,'Lz','1H')'*rho)`.

### Control flow inferred from the code

- Line 78: `for` loop over `n=1:nsteps`.
- Line 89: `for` loop over `n=2:nsteps`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=9.403`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'E','1H'}`.
- Lines 22: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.00614 2.00194 2.00988]`.
- Lines 23: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=pi*[253.6 105.1 123.8]/180`.
- Lines 24: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[0.00 0.00 0.00]`.
- Lines 25: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[0.00 0.00 0.00]`.
- Lines 26-27: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00], [0.00 0.00 3.00]}`.
- Lines 30: computes `inter.relaxation` using `inter.relaxation={'weizmann'}`.
- Lines 31: computes `inter.weiz_r1e` using `inter.weiz_r1e=1/0.3e-3`.
- Lines 32: computes `inter.weiz_r1n` using `inter.weiz_r1n=1/4.0`.
- Lines 33: computes `inter.weiz_r2e` using `inter.weiz_r2e=1/1.0e-6`.
- Lines 34: computes `inter.weiz_r2n` using `inter.weiz_r2n=1/0.2e-3`.
- Lines 35: computes `inter.weiz_r1d` using `inter.weiz_r1d=zeros(2,2)`.
- Lines 36: computes `inter.weiz_r2d` using `inter.weiz_r2d=zeros(2,2)`.
- Lines 37: computes `inter.temperature` using `inter.temperature=100`.
- Lines 38: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 39: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 42: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.

## Implementation structure

- A MAS DNP simulation performed as described in Fred Mentink-
- Vigier's paper (Spinach rotation conventions are different):
- Steady state rotor period simulation for a single crystal,
- computed using Newton-Raphson steady state solver.
- Calculation time: seconds
- Magnet field
- Spin specification
- Interactions
- Relaxation parameters
- Basis set
- Spinach housekeeping
- Stack generation parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rotor_stack()`, `operator()`, `relaxation()`, `speye()`, `propagator()`, `steady()`, `rho()`, `step()`, `kfigure()`, `trajan()`, `equilibrium()`, `hamiltonian()`, `assume()`, `state()`.
