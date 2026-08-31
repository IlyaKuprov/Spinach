# examples/nqr/pure_nqr_iodine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nqr/pure_nqr_iodine.m`
- Signature: `pure_nqr_iodine()`
- Total lines: 53

## Purpose

Powder NQR spectrum of a system with a single 127I nucleus. Calculation time: seconds

## Physical / mathematical content

- NQR examples. The Hamiltonian is dominated by quadrupolar interaction with little or no Zeeman field, so transition frequencies reflect electric field gradients and asymmetry parameters.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: System specification; implemented by `sys.magnet=0`.
- Lines 14-15: Formalism and basis; implemented by `bas.formalism='sphten-liouv'`.
- Lines 18-19: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 25-26: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 29-30: Experiment parameters; implemented by `parameters.spins={'127I'}`.
- Lines 40-41: Simulation; implemented by `fid=powder(spin_system,@hp_acquire,parameters,'labframe')`.
- Lines 43-44: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 46-47: Fourier transform; implemented by `spec=imag(fftshift(fft(fid)))`.
- Lines 49-50: Plotting; implemented by `kfigure(); plot_1d(spin_system,spec,parameters)`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=0`.
- Lines 11: computes `sys.isotopes` using `sys.isotopes={'127I'}`.
- Lines 12: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(560e6,0.01,5/2,[0 0 0])`.
- Lines 15: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 16: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 19: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 20: computes `inter.damp_rate` using `inter.damp_rate=1e5`.
- Lines 21: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 22: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 23: computes `inter.temperature` using `inter.temperature=298`.
- Lines 26: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 30: computes `parameters.spins` using `parameters.spins={'127I'}`.
- Lines 31: computes `parameters.needs` using `parameters.needs={'aniso_eq'}`.
- Lines 32: computes `parameters.sweep` using `parameters.sweep=5e8`.
- Lines 33: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 34: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 35: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','127I')`.
- Lines 36: computes `parameters.pulse_op` using `parameters.pulse_op=operator(spin_system,'Lx','127I')`.

## Implementation structure

- Powder NQR spectrum of a system with a single 127I nucleus.
- Calculation time: seconds
- System specification
- Formalism and basis
- Relaxation theory
- Spinach housekeeping
- Experiment parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `operator()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
