# examples/fundamentals/state_spaces_4.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_spaces_4.m`
- Signature: `state_spaces_4()`
- Total lines: 56

## Purpose

Trajectory analysis for a MAS simulation of isotopically labelled glycine powder, starting from L+ on protons. Calculation time: hours, faster on a GPU.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-13: Spin system properties (PCM DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/glycine.log'), {{'H','1H'},{'C','13C'},{'N','15N'}}, [31.5 182.1 264.5],[])`.
- Lines 14-15: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 17-18: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 22-23: Force Krylov propagation; implemented by `sys.tols.krylov_tol=1000`.
- Lines 25-29: This needs a GPU sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Experiment setup; implemented by `parameters.rate=2000`.
- Lines 45-46: Get the trajectory; implemented by `traj=singlerot(spin_system,@traject,parameters,'nmr')`.
- Lines 48-49: Average over rotor phase; implemented by `traj=fpl2rho(traj,2*parameters.max_rank+1)`.
- Lines 51-52: Trajectory analysis; implemented by `kfigure(); trajan(spin_system,traj,'correlation_order')`.

### Key state/data transformations

- Lines 11-13: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/glycine.log'), {{'H','1H'},{'C','13C'},{'N','15N'}}, [31.5 182.1 264.5],[])`.
- Lines 15: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 18: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 19: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 20: computes `bas.longitudinals` using `bas.longitudinals={'15N','13C'}`.
- Lines 23: computes `sys.tols.krylov_tol` using `sys.tols.krylov_tol=1000`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `parameters.rate` using `parameters.rate=2000`.
- Lines 34: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 35: computes `parameters.max_rank` using `parameters.max_rank=17`.
- Lines 36: computes `parameters.sweep` using `parameters.sweep=1e5`.
- Lines 37: computes `parameters.npoints` using `parameters.npoints=64`.
- Lines 38: computes `parameters.offset` using `parameters.offset=15000`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'13C'}`.
- Lines 40: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 41: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 42: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 43: computes `parameters.verbose` using `parameters.verbose=1`.

## Implementation structure

- Trajectory analysis for a MAS simulation of isotopically labelled
- glycine powder, starting from L+ on protons.
- Calculation time: hours, faster on a GPU.
- Spin system properties (PCM DFT calculation)
- Magnet field
- Basis set
- Force Krylov propagation
- This needs a GPU
- sys.enable={'gpu'};
- Spinach housekeeping
- Experiment setup
- Get the trajectory

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `singlerot()`, `fpl2rho()`, `kfigure()`, `trajan()`, `ylim()`, `set()`.
