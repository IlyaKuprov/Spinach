# examples/fundamentals/state_spaces_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_spaces_3.m`
- Signature: `state_spaces_3()`
- Total lines: 64

## Purpose

Transverse magnetisation dynamics in a pulse-acquire experiment on a fatty acid. This example looks at how the magnetisation drifts around the state space under the influence of strong J-coupling in the absence of relaxation. Calculation time: minutes.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Read spin system properties; implemented by `[sys,inter]=fatty_acid(15)`.
- Lines 15-16: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 18-19: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Algorithmic options; implemented by `sys.enable={'greedy','prop_cache'}`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 36-37: Assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 39-40: Operators; implemented by `H=hamiltonian(spin_system)`.
- Lines 43-44: Trajectory generation; implemented by `traj=evolution(spin_system,H,[],parameters.rho0,4e-5,50,'trajectory')`.
- Lines 46-47: Number of repetition; implemented by `N_rep=8`.
- Lines 49-50: CPMG loop; implemented by `for n=1:N_rep`.
- Lines 52-53: Pulse; implemented by `traj(:,end)=step(spin_system,Lx,traj(:,end),pi)`.
- Lines 55-56: Trajectory generation; implemented by `traj=[traj evolution(spin_system,H,[],traj(:,end),4e-5,100,'trajectory')]`.
- Lines 60-61: Trajectory analysis; implemented by `kfigure(); trajan(spin_system,traj,'correlation_order')`.

### Control flow inferred from the code

- Line 50: `for` loop over `n=1:N_rep`.

### Key state/data transformations

- Lines 13: computes `[sys,inter]` using `[sys,inter]=fatty_acid(15)`.
- Lines 16: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 19: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 20: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 21: computes `bas.space_level` using `bas.space_level=1`.
- Lines 22: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 25: computes `sys.enable` using `sys.enable={'greedy','prop_cache'}`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 33: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lx','1H')`.
- Lines 34: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 40: computes `H` using `H=hamiltonian(spin_system)`.
- Lines 41: computes `Lx` using `Lx=operator(spin_system,'Lx','1H')`.
- Lines 44: computes `traj` using `traj=evolution(spin_system,H,[],parameters.rho0,4e-5,50,'trajectory')`.
- Lines 47: computes `N_rep` using `N_rep=8`.
- Lines 53: computes `traj(:,end)` using `traj(:,end)=step(spin_system,Lx,traj(:,end),pi)`.

## Implementation structure

- Transverse magnetisation dynamics in a pulse-acquire experiment on
- a fatty acid. This example looks at how the magnetisation drifts
- around the state space under the influence of strong J-coupling in
- the absence of relaxation.
- Calculation time: minutes.
- Read spin system properties
- Magnet field
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Assumptions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `fatty_acid()`, `create()`, `basis()`, `state()`, `assume()`, `hamiltonian()`, `operator()`, `evolution()`, `traj()`, `step()`, `kfigure()`, `trajan()`.
