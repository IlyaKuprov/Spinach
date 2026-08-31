# examples/relaxation_theory/inv_rec_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/inv_rec_1.m`
- Signature: `inv_rec_1()`
- Total lines: 63

## Purpose

A simple inversion-recovery experiment; longitudinal magnetisation is monitored as a function of time. Calculation time: seconds.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Spin system; implemented by `sys.magnet=14.1`.
- Lines 14-15: Zeeman interactions; implemented by `inter.zeeman.scalar={1.5}`.
- Lines 17-18: Complete basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 21-22: Relaxation theory; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 29-30: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Isotropic thermal equilibrium; implemented by `rho=equilibrium(spin_system)`.
- Lines 36-37: Detection state; implemented by `coil=state(spin_system,'Lz','1H')`.
- Lines 39-40: Static Liouvillian superoperator; implemented by `L=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 42-43: Pulse operator; implemented by `Lx=operator(spin_system,'Lx','1H')`.
- Lines 45-46: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 48-49: Inversion pulse; implemented by `rho=step(spin_system,Lx,rho,pi)`.
- Lines 51-52: Evolution (1000 steps, 1.0 ms each); implemented by `answer=evolution(spin_system,L+1i*R,coil,rho,1e-3,1000,'observable')`.
- Lines 54-55: Plotting; implemented by `kfigure(); scale_figure([1.00 0.65])`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.5}`.
- Lines 18: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 19: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 22: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 23: computes `inter.r1_rates` using `inter.r1_rates={5.0}`.
- Lines 24: computes `inter.r2_rates` using `inter.r2_rates={5.0}`.
- Lines 25: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 27: computes `inter.temperature` using `inter.temperature=298`.
- Lines 30: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `rho` using `rho=equilibrium(spin_system)`.
- Lines 37: computes `coil` using `coil=state(spin_system,'Lz','1H')`.
- Lines 40: computes `L` using `L=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 43: computes `Lx` using `Lx=operator(spin_system,'Lx','1H')`.
- Lines 46: computes `R` using `R=relaxation(spin_system)`.
- Lines 52: computes `answer` using `answer=evolution(spin_system,L+1i*R,coil,rho,1e-3,1000,'observable')`.

## Implementation structure

- A simple inversion-recovery experiment; longitudinal
- magnetisation is monitored as a function of time.
- Calculation time: seconds.
- Spin system
- Zeeman interactions
- Complete basis set
- Relaxation theory
- Spinach housekeeping
- Isotropic thermal equilibrium
- Detection state
- Static Liouvillian superoperator
- Pulse operator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `equilibrium()`, `state()`, `hamiltonian()`, `assume()`, `operator()`, `relaxation()`, `step()`, `evolution()`, `kfigure()`, `scale_figure()`, `kylabel()`, `kxlabel()`.
