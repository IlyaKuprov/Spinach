# examples/relaxation_theory/from_md/sucrose_three_spins.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/from_md/sucrose_three_spins.m`
- Signature: `sucrose_three_spins()`
- Total lines: 112

## Purpose

One of the calculations reported in the JMR paper with Jim Prestegard: a three-spin subsystem from the glucose ring of sucrose -an illustra- tion of incorrect viscosity of TIP3P water. The experimental value of tau_c for sucrose is around 90 ps (and this is correctly reproduced by OPC and TIP5P water), but TIP3P only agrees with Redfield theory when tau_c is set to 37 ps in the latter. Here, numerical relaxation supe

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 23-24: Three-spin system; implemented by `sys.magnet=14.1`.
- Lines 27-28: Complete basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Analytical relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 38-39: Run Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Set laboratory frame assumptions; implemented by `spin_system=assume(spin_system,'labframe')`.
- Lines 45-46: Get coherent Hamiltonian; implemented by `H0=hamiltonian(spin_system)`.
- Lines 48-49: Load MD trajectory, dims: (XYZ, spins, time); implemented by `load('suc_three_spin_traj.mat','traj','dt')`.
- Lines 51-52: 50k frames is enough here; implemented by `traj=traj(:,:,1:50000); ref_frame=1`.
- Lines 54-55: Compute Hamiltonian for each trajectory frame; implemented by `report(spin_system,'computing MD frame Hamiltonians ')`.
- Lines 59-60: Localise system object; implemented by `spinsys_loc=spin_system`.
- Lines 62-63: Hush up Spinach output; implemented by `spinsys_loc.sys.output='hush'`.
- Lines 65-66: Pull a trajectory frame; implemented by `traj_slice=traj(:,:,n)`.
- Lines 68-69: Set current coordinates; implemented by `spinsys_loc.inter.coordinates={double(traj_slice(:,1)')`.
- Lines 74-75: Get the anisotropic part of the Hamiltonian; implemented by `[~,Q]=hamiltonian(spinsys_loc); H1{n}=orientation(Q,[0,0,0])`.
- Lines 79-80: Get relaxation superoperator from the MD trajectory; implemented by `[R_gce,dR_gce]=ngce(spin_system,H0,H1,dt,inter.tau_c{1},0)`.
- Lines 82-83: Load reference frame coordinates; implemented by `spin_system.inter.coordinates={double(traj(:,1,ref_frame)')`.
- Lines 88-89: Get analytical Redfield matrix; implemented by `R_red=relaxation(spin_system)`.
- Lines 91-92: Get relaxation rates; implemented by `gce_rates_stdm=diag(dR_gce)`.

### Control flow inferred from the code

- Line 57: `parfor` loop over `n=1:size(traj,3)`.

### Key state/data transformations

- Lines 24: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 25: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H'}`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 33: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 34: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 35: computes `inter.tau_c` using `inter.tau_c={37e-12}`.
- Lines 36: computes `inter.temperature` using `inter.temperature=298`.
- Lines 39: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 46: computes `H0` using `H0=hamiltonian(spin_system)`.
- Lines 52: computes `traj` using `traj=traj(:,:,1:50000); ref_frame=1`.
- Lines 56: computes `H1` using `H1=cell(1,size(traj,3))`.
- Lines 60: computes `spinsys_loc` using `spinsys_loc=spin_system`.
- Lines 63: computes `spinsys_loc.sys.output` using `spinsys_loc.sys.output='hush'`.
- Lines 66: computes `traj_slice` using `traj_slice=traj(:,:,n)`.
- Lines 69: computes `spinsys_loc.inter.coordinates` using `spinsys_loc.inter.coordinates={double(traj_slice(:,1)')`.
- Lines 72: computes `spinsys_loc.inter.pbc` using `spinsys_loc.inter.pbc={}; spinsys_loc=dipolar(spinsys_loc)`.

## Implementation structure

- One of the calculations reported in the JMR paper with Jim Prestegard:
- a three-spin subsystem from the glucose ring of sucrose -an illustra-
- tion of incorrect viscosity of TIP3P water.
- The experimental value of tau_c for sucrose is around 90 ps (and this
- is correctly reproduced by OPC and TIP5P water), but TIP3P only agrees
- with Redfield theory when tau_c is set to 37 ps in the latter.
- Here, numerical relaxation superoperator computed from a long MD tra-
- jectory is compared with the analytical one computed using the isotro-
- pic rotational diffusion approximation.
- Calculation time: minutes, with most of the time spent
- computing MD frame Hamiltonians
- Three-spin system

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `hamiltonian()`, `load()`, `traj()`, `report()`, `double()`, `traj_slice()`, `dipolar()`, `orientation()`, `ngce()`, `relaxation()`, `kfigure()`, `errorbar()`, `kxlabel()`.
