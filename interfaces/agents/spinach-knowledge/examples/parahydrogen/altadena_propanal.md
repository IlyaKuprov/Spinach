# examples/parahydrogen/altadena_propanal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/parahydrogen/altadena_propanal.m`
- Signature: `altadena_propanal()`
- Total lines: 68

## Purpose

ALTADENA experiment simulation for the parahydrogenation of acrolein into propanal. Simple model of the ALTADENA effect is used: perfect- ly adiabatic transfer is assumed and the isotropic mixing in low fi- eld is ignored completely. Note the small flip angle. Calculation time: seconds

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Spin system; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 16-17: Magnetic field; implemented by `sys.magnet=7.05`.
- Lines 19-20: Chemical shifts; implemented by `inter.zeeman.scalar={1.11 1.11 1.11 2.46 2.46 9.79}`.
- Lines 22-23: Scalar couplings; implemented by `inter.coupling.scalar=cell(6)`.
- Lines 29-30: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 55-56: Simulation; implemented by `fid=liquid(spin_system,@hp_acquire,parameters,'nmr')`.
- Lines 58-59: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 61-62: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 64-65: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 17: computes `sys.magnet` using `sys.magnet=7.05`.
- Lines 20: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.11 1.11 1.11 2.46 2.46 9.79}`.
- Lines 23: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(6)`.
- Lines 24: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=7.3; inter.coupling.scalar{2,4}=7.3`.
- Lines 25: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=7.3; inter.coupling.scalar{1,5}=7.3`.
- Lines 26: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}=7.3; inter.coupling.scalar{3,5}=7.3`.
- Lines 27: computes `inter.coupling.scalar{4,6}` using `inter.coupling.scalar{4,6}=1.4; inter.coupling.scalar{5,6}=1.4`.
- Lines 30: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 31: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `bas.sym_group` using `bas.sym_group={'S3','S2'}`.
- Lines 33: computes `bas.sym_spins` using `bas.sym_spins={[1 2 3],[4 5]}`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 41-43: computes `parameters.rho0` using `parameters.rho0=1.0*state(spin_system,{'Lz','Lz'},{1,4})- 0.5*state(spin_system,{'Lz'},{1})+ 0.5*state(spin_system,{'Lz'},{4})`.
- Lines 44: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 45: computes `parameters.pulse_op` using `parameters.pulse_op=operator(spin_system,'Ly','1H')`.
- Lines 46: computes `parameters.pulse_angle` using `parameters.pulse_angle=pi/100`.

## Implementation structure

- ALTADENA experiment simulation for the parahydrogenation of acrolein
- into propanal. Simple model of the ALTADENA effect is used: perfect-
- ly adiabatic transfer is assumed and the isotropic mixing in low fi-
- eld is ignored completely. Note the small flip angle.
- Calculation time: seconds
- Spin system
- Magnetic field
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
