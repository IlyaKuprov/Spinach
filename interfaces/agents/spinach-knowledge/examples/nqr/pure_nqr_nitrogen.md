# examples/nqr/pure_nqr_nitrogen.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nqr/pure_nqr_nitrogen.m`
- Signature: `pure_nqr_nitrogen()`
- Total lines: 54

## Purpose

Powder NQR spectrum of a system with a single 14N nucleus. Calculation time: seconds

## Physical / mathematical content

- NQR examples. The Hamiltonian is dominated by quadrupolar interaction with little or no Zeeman field, so transition frequencies reflect electric field gradients and asymmetry parameters.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: System specification; implemented by `sys.magnet=0`.
- Lines 15-16: Formalism and basis; implemented by `bas.formalism='sphten-liouv'`.
- Lines 19-20: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 26-27: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Experiment parameters; implemented by `parameters.spins={'14N'}`.
- Lines 41-42: Simulation; implemented by `fid=powder(spin_system,@hp_acquire,parameters,'labframe')`.
- Lines 44-45: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 47-48: Fourier transform; implemented by `spec=imag(fftshift(fft(fid)))`.
- Lines 50-51: Plotting; implemented by `kfigure(); plot_1d(spin_system,spec,parameters)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=0`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'14N'}`.
- Lines 13: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.53,1,[0 0 0])`.
- Lines 16: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 17: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 20: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 21: computes `inter.damp_rate` using `inter.damp_rate=1e5`.
- Lines 22: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 23: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 24: computes `inter.temperature` using `inter.temperature=298`.
- Lines 27: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `parameters.spins` using `parameters.spins={'14N'}`.
- Lines 32: computes `parameters.needs` using `parameters.needs={'aniso_eq'}`.
- Lines 33: computes `parameters.sweep` using `parameters.sweep=5e6`.
- Lines 34: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 35: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 36: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','14N')`.
- Lines 37: computes `parameters.pulse_op` using `parameters.pulse_op=operator(spin_system,'Lx','14N')`.

## Implementation structure

- Powder NQR spectrum of a system with a single 14N nucleus.
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
