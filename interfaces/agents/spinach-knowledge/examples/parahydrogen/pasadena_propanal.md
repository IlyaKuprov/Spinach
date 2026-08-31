# examples/parahydrogen/pasadena_propanal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/parahydrogen/pasadena_propanal.m`
- Signature: `pasadena_propanal()`
- Total lines: 64

## Purpose

PASADENA experiment simulation for the parahydrogenation of acrolein into propanal. Calculation time: seconds

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Spin system; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 14-15: Magnetic field; implemented by `sys.magnet=7.05`.
- Lines 17-18: Chemical shifts; implemented by `inter.zeeman.scalar={1.11 1.11 1.11 2.46 2.46 9.79}`.
- Lines 20-21: Scalar couplings; implemented by `inter.coupling.scalar=cell(6)`.
- Lines 27-28: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 33-34: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 37-38: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 51-52: Simulation; implemented by `fid=liquid(spin_system,@hp_acquire,parameters,'nmr')`.
- Lines 54-55: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'gauss',10}})`.
- Lines 57-58: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 60-61: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 15: computes `sys.magnet` using `sys.magnet=7.05`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.11 1.11 1.11 2.46 2.46 9.79}`.
- Lines 21: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(6)`.
- Lines 22: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=7.3; inter.coupling.scalar{2,4}=7.3`.
- Lines 23: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=7.3; inter.coupling.scalar{1,5}=7.3`.
- Lines 24: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}=7.3; inter.coupling.scalar{3,5}=7.3`.
- Lines 25: computes `inter.coupling.scalar{4,6}` using `inter.coupling.scalar{4,6}=1.4; inter.coupling.scalar{5,6}=1.4`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `bas.sym_group` using `bas.sym_group={'S3','S2'}`.
- Lines 31: computes `bas.sym_spins` using `bas.sym_spins={[1 2 3],[4 5]}`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 38: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 39: computes `parameters.rho0` using `parameters.rho0=state(spin_system,{'Lz','Lz'},{1,4})`.
- Lines 40: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 41: computes `parameters.pulse_op` using `parameters.pulse_op=operator(spin_system,'Ly','1H')`.
- Lines 42: computes `parameters.pulse_angle` using `parameters.pulse_angle=pi/4`.

## Implementation structure

- PASADENA experiment simulation for the parahydrogenation of acrolein
- into propanal.
- Calculation time: seconds
- Spin system
- Magnetic field
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
