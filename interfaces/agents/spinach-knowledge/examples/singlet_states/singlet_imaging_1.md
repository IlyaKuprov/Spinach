# examples/singlet_states/singlet_imaging_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/singlet_imaging_1.m`
- Signature: `singlet_imaging_1()`
- Total lines: 160

## Purpose

Singlet imaging in a system with one-dimensional diffusion and flow. Calculation time: minutes

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file also defines local helper function(s): `tube_flow()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Spin system and interactions; implemented by `sys.magnet=9.4`.
- Lines 21-22: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 27-28: Relaxation superoperator accuracy; implemented by `sys.tols.rlx_integration=1e-5`.
- Lines 31-32: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 34-35: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 38-39: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Sample geometry; implemented by `parameters.dims=[0.10 0.015]`.
- Lines 47-48: Sequence parameters; implemented by `parameters.spins={'13C'}`.
- Lines 51-52: Assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 54-55: Relaxation phantom; implemented by `parameters.rlx_ph={ones(parameters.npts)}`.
- Lines 58-59: NMR sample tube phantom; implemented by `tube=zeros(parameters.npts); tube(:,6:10)=1`.
- Lines 61-62: Initial and detection state phantoms; implemented by `parameters.rho0_ph={tube}`.
- Lines 67-68: Diffusion and flow; implemented by `parameters.u=-6e-2*ones(parameters.npts)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'13C','13C'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.03,-0.03}`.
- Lines 16: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2)`.
- Lines 17: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=55`.
- Lines 18: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 22: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 23: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 24: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 25: computes `inter.tau_c` using `inter.tau_c={1e-9}`.
- Lines 28: computes `sys.tols.rlx_integration` using `sys.tols.rlx_integration=1e-5`.
- Lines 29: computes `sys.tols.rlx_zero` using `sys.tols.rlx_zero=1e-5`.
- Lines 32: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 35: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 36: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 39: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 43: computes `parameters.dims` using `parameters.dims=[0.10 0.015]`.
- Lines 44: computes `parameters.npts` using `parameters.npts=[150 15]`.

### Local helper functions

- Line 73: `tube_flow()` — `function traj=tube_flow(spin_system,parameters,H,R,K,G,F)`. Compose Liouvillian
  - Representative operation: `L=H+F+1i*R+1i*K`.
  - Representative operation: `Lx=operator(spin_system,'Lx',parameters.spins{1})`.

## Implementation structure

- Singlet imaging in a system with one-dimensional
- diffusion and flow.
- Calculation time: minutes
- Spin system and interactions
- Relaxation theory
- Relaxation superoperator accuracy
- Algorithmic options
- Basis set
- Spinach housekeeping
- Sample geometry
- Sequence parameters
- Assumptions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `relaxation()`, `tube()`, `state()`, `tube_flow()`, `operator()`, `speye()`, `pulse_shape()`, `shaped_pulse_af()`, `step()`, `evolution()`, `imaging()`, `singlet()`, `kfigure()`.
