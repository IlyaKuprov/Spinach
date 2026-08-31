# examples/parahydrogen/case_studies/spontaneous_singlet_to_z/rlx_trajectory.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/parahydrogen/case_studies/spontaneous_singlet_to_z/rlx_trajectory.m`
- Signature: `rlx_trajectory()`
- Total lines: 75

## Purpose

Time dependence of LzSz and Lz+Sz spin orders in a para- hydrogen molecule coordinated to a nickel cage that cre- ates large chemical shift anisotropy. Done for Gloggler group, paper link coming in due course; see also the Ma- thematica worksheet.

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnet field, Tesla; implemented by `sys.magnet=18.7893`.
- Lines 14-15: Isotopes; implemented by `sys.isotopes={'1H','1H'}`.
- Lines 17-18: Coordinates (for dipole tensor); implemented by `inter.coordinates={[ 0.1399 16.6491 21.2341]`.
- Lines 21-22: Zeeman interaction tensors, traces subtracted; implemented by `inter.zeeman.matrix={[61.897647 13.533565 1.607280`.
- Lines 31-32: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 37-38: Formalism and basis; implemented by `bas.formalism='sphten-liouv'`.
- Lines 41-42: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 45-46: Assumptions; implemented by `spin_system=assume(spin_system,'nmr')`.
- Lines 48-49: Build the Liouvillian; implemented by `L=hamiltonian(spin_system)+1i*relaxation(spin_system)`.
- Lines 51-52: Initial state; implemented by `rho=singlet(spin_system,1,2)`.
- Lines 54-56: Detection states; implemented by `A=state(spin_system,{'L+','L-'},{1 2})+ state(spin_system,{'L-','L+'},{1 2})`.
- Lines 60-61: Time evolution; implemented by `traj=evolution(spin_system,L,[A/2 4*B C],rho,2e-3,1000,'multichannel')`.
- Lines 63-64: Plotting; implemented by `time_axis=linspace(0,2,1001)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=18.7893`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 18: computes `inter.coordinates` using `inter.coordinates={[ 0.1399 16.6491 21.2341]`.
- Lines 22: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={[61.897647 13.533565 1.607280`.
- Lines 32: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 33: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 34: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 35: computes `inter.tau_c` using `inter.tau_c={500e-12}`.
- Lines 38: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 39: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 42: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 49: computes `L` using `L=hamiltonian(spin_system)+1i*relaxation(spin_system)`.
- Lines 52: computes `rho` using `rho=singlet(spin_system,1,2)`.
- Lines 55-56: computes `A` using `A=state(spin_system,{'L+','L-'},{1 2})+ state(spin_system,{'L-','L+'},{1 2})`.
- Lines 57: computes `B` using `B=state(spin_system,{'Lz','Lz'},{1 2})`.
- Lines 58: computes `C` using `C=state(spin_system,'Lz','all')`.
- Lines 61: computes `traj` using `traj=evolution(spin_system,L,[A/2 4*B C],rho,2e-3,1000,'multichannel')`.
- Lines 64: computes `time_axis` using `time_axis=linspace(0,2,1001)`.

## Implementation structure

- Time dependence of LzSz and Lz+Sz spin orders in a para-
- hydrogen molecule coordinated to a nickel cage that cre-
- ates large chemical shift anisotropy. Done for Gloggler
- group, paper link coming in due course; see also the Ma-
- thematica worksheet.
- Magnet field, Tesla
- Isotopes
- Coordinates (for dipole tensor)
- Zeeman interaction tensors, traces subtracted
- Relaxation theory parameters
- Formalism and basis
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cellfun()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `relaxation()`, `singlet()`, `state()`, `evolution()`, `kfigure()`, `kxlabel()`, `kylabel()`, `klegend()`, `scale_figure()`.
