# examples/parahydrogen/pasadena_ethylbenzene.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/parahydrogen/pasadena_ethylbenzene.m`
- Signature: `pasadena_ethylbenzene()`
- Total lines: 70

## Purpose

PASADENA experiment simulation for the parahydrogenation of styrene into ethylbenzene. Set to reproduce the top trace of Fig 5 in Calculation time: seconds

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Spin system; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 15-16: Magnetic field; implemented by `sys.magnet=7.05`.
- Lines 18-20: Chemical shifts; implemented by `inter.zeeman.scalar={1.201 1.201 1.201 2.625 2.625 7.207 7.207 7.265 7.265 7.155}`.
- Lines 22-23: Scalar couplings; implemented by `inter.coupling.scalar=cell(10)`.
- Lines 31-32: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 39-40: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 43-44: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 57-58: Simulation; implemented by `fid=liquid(spin_system,@hp_acquire,parameters,'nmr')`.
- Lines 60-61: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'gauss',10}})`.
- Lines 63-64: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 66-67: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 16: computes `sys.magnet` using `sys.magnet=7.05`.
- Lines 19-20: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.201 1.201 1.201 2.625 2.625 7.207 7.207 7.265 7.265 7.155}`.
- Lines 23: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(10)`.
- Lines 24: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=7.63; inter.coupling.scalar{2,4}=7.63`.
- Lines 25: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=7.63; inter.coupling.scalar{1,5}=7.63`.
- Lines 26: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}=7.63; inter.coupling.scalar{3,5}=7.63`.
- Lines 27: computes `inter.coupling.scalar{6,8}` using `inter.coupling.scalar{6,8}=7.63; inter.coupling.scalar{6,10}=1.26`.
- Lines 28: computes `inter.coupling.scalar{7,9}` using `inter.coupling.scalar{7,9}=7.63; inter.coupling.scalar{7,10}=1.26`.
- Lines 29: computes `inter.coupling.scalar{8,10}` using `inter.coupling.scalar{8,10}=7.44; inter.coupling.scalar{9,10}=7.44`.
- Lines 32: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 33: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 34: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 35: computes `bas.space_level` using `bas.space_level=1`.
- Lines 36: computes `bas.sym_group` using `bas.sym_group={'S3','S2'}`.
- Lines 37: computes `bas.sym_spins` using `bas.sym_spins={[1 2 3],[4 5]}`.
- Lines 40: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 44: computes `parameters.spins` using `parameters.spins={'1H'}`.

## Implementation structure

- PASADENA experiment simulation for the parahydrogenation of styrene
- into ethylbenzene. Set to reproduce the top trace of Fig 5 in
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
