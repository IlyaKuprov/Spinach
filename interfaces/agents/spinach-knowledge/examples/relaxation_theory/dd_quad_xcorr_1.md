# examples/relaxation_theory/dd_quad_xcorr_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/dd_quad_xcorr_1.m`
- Signature: `dd_quad_xcorr_1()`
- Total lines: 60

## Purpose

Complete Bloch-Redfield-Wangsness relaxation superoperator in a system with a quadrupolar coupling and a dipole coupling. Spinach relaxation theory module automatically accounts for all cross-correlations (dipole- quadrupole cross-correlation is present in this case). Dipolar couplings are computed from Cartesian coordinates of the two spins. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: System specification; implemented by `sys.magnet=14.1`.
- Lines 22-23: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 28-29: Basis specification; implemented by `bas.formalism='sphten-liouv'`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 47-48: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 50-51: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 53-54: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 56-57: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'1H','14N'}`.
- Lines 16: computes `inter.coupling.scalar` using `inter.coupling.scalar={0 50; 50 0}`.
- Lines 17: computes `inter.coupling.eigs{2,2}` using `inter.coupling.eigs{2,2}=[1e4 1e4 -2e4]`.
- Lines 18: computes `inter.coupling.euler{2,2}` using `inter.coupling.euler{2,2}=[0 0 0]`.
- Lines 19: computes `inter.coordinates` using `inter.coordinates={[0 0 0]`.
- Lines 23: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 24: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 25: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 26: computes `inter.tau_c` using `inter.tau_c={1e-9}`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 38: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 39: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 40: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 41: computes `parameters.offset` using `parameters.offset=0`.

## Implementation structure

- Complete Bloch-Redfield-Wangsness relaxation superoperator in a system
- with a quadrupolar coupling and a dipole coupling. Spinach relaxation
- theory module automatically accounts for all cross-correlations (dipole-
- quadrupole cross-correlation is present in this case). Dipolar couplings
- are computed from Cartesian coordinates of the two spins.
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Basis specification
- Spinach housekeeping
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
