# examples/spin_chemistry/cidnp_flash_acquire.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/cidnp_flash_acquire.m`
- Signature: `cidnp_flash_acquire()`
- Total lines: 100

## Purpose

A model of the CIDNP magnetisation pumping process described in IK's paper: The system is pumped for 0.5 seconds, and then allowed to relax. Calculation time: seconds. Miguel Mompean Ilya Kuprov

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 18-19: Isotopes; implemented by `sys.isotopes={'1H','19F'}`.
- Lines 21-22: Chemical shifts; implemented by `inter.zeeman.scalar={0.0, 0.0}`.
- Lines 24-25: Chemical shift anisotropies (DFT); implemented by `inter.zeeman.eigs{1}=[0 0 0]`.
- Lines 30-32: Coordinates (DFT); implemented by `inter.coordinates={[0.00 0.00 0.00], [0.00 2.60 0.00]}`.
- Lines 34-35: J-coupling (expt); implemented by `inter.coupling.scalar=cell(2,2)`.
- Lines 38-39: Relaxation theories; implemented by `inter.relaxation={'redfield'}`.
- Lines 44-45: Formalism and basis; implemented by `bas.formalism='sphten-liouv'`.
- Lines 48-49: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 52-53: Get the relevant states; implemented by `Hz=state(spin_system,{'Lz'},{1})`.
- Lines 59-60: Get the Hamiltonian; implemented by `H=hamiltonian(assume(spin_system,'nmr'))`.
- Lines 62-63: Get relaxation matrix; implemented by `R=relaxation(spin_system)`.
- Lines 65-66: Thermalise R to Hz+Fz; implemented by `R(1,1)=1; R(:,1)=-R*(Hz+Fz)`.
- Lines 68-69: Add pumping terms; implemented by `R_light=magpump(spin_system,R,Hz,1.3)`.
- Lines 72-73: Start at equilibrium; implemented by `rho0=unit_state(spin_system)+Hz+Fz`.
- Lines 75-76: Run the evolution with illumination, 50 steps of 0.01 seconds; implemented by `traj_a=evolution(spin_system,H+1i*R_light,[],rho0,0.01,50,'trajectory')`.
- Lines 78-79: Run the evolution with illumination, 500 steps of 0.01 seconds; implemented by `traj_b=evolution(spin_system,H+1i*R,[],traj_a(:,end),0.01,500,'trajectory')`.
- Lines 81-82: Concatenate trajectories; implemented by `traj=[traj_a traj_b(:,2:end)]`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'1H','19F'}`.
- Lines 22: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0, 0.0}`.
- Lines 25: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[0 0 0]`.
- Lines 26: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[0 0 0]`.
- Lines 27: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[-47 -16 63]`.
- Lines 28: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[0 0 0]`.
- Lines 31-32: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00], [0.00 2.60 0.00]}`.
- Lines 35: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2,2)`.
- Lines 36: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=50`.
- Lines 39: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 40: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 41: computes `inter.tau_c` using `inter.tau_c={110e-12}`.
- Lines 42: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 45: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 46: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 49: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 53: computes `Hz` using `Hz=state(spin_system,{'Lz'},{1})`.

## Implementation structure

- A model of the CIDNP magnetisation pumping process described in
- IK's paper:
- The system is pumped for 0.5 seconds, and then allowed to relax.
- Calculation time: seconds.
- Miguel Mompean
- Ilya Kuprov
- Magnet field
- Isotopes
- Chemical shifts
- Chemical shift anisotropies (DFT)
- Coordinates (DFT)
- J-coupling (expt)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `hamiltonian()`, `assume()`, `relaxation()`, `magpump()`, `unit_state()`, `evolution()`, `traj_a()`, `traj_b()`, `kfigure()`, `scale_figure()`, `subplot()`, `ktitle()`, `kxlabel()`.
