# examples/relaxation_theory/sat_rec_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/sat_rec_1.m`
- Signature: `sat_rec_1()`
- Total lines: 61

## Purpose

A simple saturation-recovery experiment. Calculation time: seconds.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Spin system; implemented by `sys.magnet=14.1`.
- Lines 13-14: Zeeman interactions; implemented by `inter.zeeman.scalar={1.5}`.
- Lines 16-17: Complete basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 20-21: Relaxation theory; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Initial state -unit; implemented by `rho=state(spin_system,'E','1H')`.
- Lines 35-36: Detection state; implemented by `coil=state(spin_system,'Lz','1H')`.
- Lines 38-39: Static Hamiltonian superoperator; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 41-42: Pulse operator; implemented by `Lx=operator(spin_system,'Lx','1H')`.
- Lines 44-45: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 47-48: Inversion pulse; implemented by `rho=step(spin_system,Lx,rho,pi)`.
- Lines 50-51: Evolution (1000 steps, 1.0 ms each); implemented by `answer=evolution(spin_system,H+1i*R,coil,rho,1e-3,1000,'observable')`.
- Lines 53-54: Plotting; implemented by `kfigure(); x_axis=linspace(0,1,1001)`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 11: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 14: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.5}`.
- Lines 17: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 18: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 21: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 22: computes `inter.r1_rates` using `inter.r1_rates={5.0}`.
- Lines 23: computes `inter.r2_rates` using `inter.r2_rates={5.0}`.
- Lines 24: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 25: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 26: computes `inter.temperature` using `inter.temperature=298`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `rho` using `rho=state(spin_system,'E','1H')`.
- Lines 36: computes `coil` using `coil=state(spin_system,'Lz','1H')`.
- Lines 39: computes `H` using `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 42: computes `Lx` using `Lx=operator(spin_system,'Lx','1H')`.
- Lines 45: computes `R` using `R=relaxation(spin_system)`.
- Lines 51: computes `answer` using `answer=evolution(spin_system,H+1i*R,coil,rho,1e-3,1000,'observable')`.

## Implementation structure

- A simple saturation-recovery experiment.
- Calculation time: seconds.
- Spin system
- Zeeman interactions
- Complete basis set
- Relaxation theory
- Spinach housekeeping
- Initial state -unit
- Detection state
- Static Hamiltonian superoperator
- Pulse operator
- Relaxation superoperator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `hamiltonian()`, `assume()`, `operator()`, `relaxation()`, `step()`, `evolution()`, `kfigure()`, `kylabel()`, `kxlabel()`.
