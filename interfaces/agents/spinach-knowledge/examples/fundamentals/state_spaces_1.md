# examples/fundamentals/state_spaces_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_spaces_1.m`
- Signature: `state_spaces_1()`
- Total lines: 60

## Purpose

Correlation order dynamics in a pulse-acquire experiment on strychnine. Set to reproduce Figure 4 from our state space restriction accuracy analysis paper: Run time: hours (much faster on a Tesla A100 GPU)

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Read spin system properties; implemented by `[sys,inter]=strychnine({'1H'})`.
- Lines 16-17: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 26-27: Proximity cut-off; implemented by `sys.tols.prox_cutoff=4.0`.
- Lines 29-30: Algorithmic options; implemented by `sys.disable={'trajlevel'}`.
- Lines 33-34: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 39-40: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 43-44: Initial condition; implemented by `rho=state(spin_system,'L+','1H')`.
- Lines 46-47: Assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 49-50: Liouvillian; implemented by `L=hamiltonian(spin_system)+1i*relaxation(spin_system)`.
- Lines 52-53: Trajectory generation; implemented by `traj=evolution(spin_system,L,[],rho,1e-3,1000,'trajectory')`.
- Lines 55-56: Trajectory analysis; implemented by `kfigure(); trajan(spin_system,traj,'correlation_order')`.

### Key state/data transformations

- Lines 14: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H'})`.
- Lines 17: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 22: computes `bas.level` using `bas.level=7; bas.space_level=1`.
- Lines 23: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 24: computes `bas.projections` using `bas.projections=1`.
- Lines 27: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 30: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 31: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 34: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 35: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 36: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 37: computes `inter.tau_c` using `inter.tau_c={200e-12}`.
- Lines 40: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 44: computes `rho` using `rho=state(spin_system,'L+','1H')`.
- Lines 50: computes `L` using `L=hamiltonian(spin_system)+1i*relaxation(spin_system)`.
- Lines 53: computes `traj` using `traj=evolution(spin_system,L,[],rho,1e-3,1000,'trajectory')`.

## Implementation structure

- Correlation order dynamics in a pulse-acquire experiment on
- strychnine. Set to reproduce Figure 4 from our state space
- restriction accuracy analysis paper:
- Run time: hours (much faster on a Tesla A100 GPU)
- Read spin system properties
- Magnet field
- Basis set
- Proximity cut-off
- Algorithmic options
- Relaxation theory parameters
- Spinach housekeeping
- Initial condition

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `state()`, `assume()`, `hamiltonian()`, `relaxation()`, `evolution()`, `kfigure()`, `trajan()`, `set()`.
