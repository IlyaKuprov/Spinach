# examples/esr_liq_pulsed/pulse_acquire_biaryl.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/pulse_acquire_biaryl.m`
- Signature: `pulse_acquire_biaryl()`
- Total lines: 85

## Purpose

A time-domain pulse-acquire version of the EasySpin biaryl test file, with acknowledgements to Stefan Stoll. The Spinach simulation is run using explicit time propagation in Liou- ville space. Full symmetry treatment using the S2xS2xS2xS2xS2xS2 group direct product is performed. Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 15-16: Magnet induction; implemented by `sys.magnet=0.33`.
- Lines 18-20: Isotope list; implemented by `sys.isotopes={'E','14N','1H','1H','1H','1H','1H', '14N','1H','1H','1H','1H','1H'}`.
- Lines 22-23: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Zeeman interactions; implemented by `inter.zeeman.scalar=cell(13,1)`.
- Lines 34-35: Spin-spin couplings; implemented by `inter.coupling.scalar=cell(13,13)`.
- Lines 49-50: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 55-56: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 59-60: Experiment parameters; implemented by `parameters.spins={'E'}`.
- Lines 72-73: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'esr')`.
- Lines 75-76: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 78-79: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 81-82: Plot the spectrum; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 16: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 19-20: computes `sys.isotopes` using `sys.isotopes={'E','14N','1H','1H','1H','1H','1H', '14N','1H','1H','1H','1H','1H'}`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `bas.longitudinals` using `bas.longitudinals={'1H','14N'}`.
- Lines 26: computes `bas.projections` using `bas.projections=+1`.
- Lines 27: computes `bas.sym_group` using `bas.sym_group={'S2','S2','S2','S2','S2','S2'}`.
- Lines 28: computes `bas.sym_spins` using `bas.sym_spins={[2 8],[3 9],[4 10],[5 11],[6 12],[7 13]}`.
- Lines 31: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=cell(13,1)`.
- Lines 32: computes `inter.zeeman.scalar{1}` using `inter.zeeman.scalar{1}=2.00316`.
- Lines 35: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(13,13)`.
- Lines 36: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2} = 12.16e6`.
- Lines 37: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3} = -6.7e6`.
- Lines 38: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4} = -1.82e6`.
- Lines 39: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5} = -7.88e6`.
- Lines 40: computes `inter.coupling.scalar{1,6}` using `inter.coupling.scalar{1,6} = -0.64e6`.
- Lines 41: computes `inter.coupling.scalar{1,7}` using `inter.coupling.scalar{1,7} = 67.93e6`.
- Lines 42: computes `inter.coupling.scalar{1,8}` using `inter.coupling.scalar{1,8} = 12.16e6`.

## Implementation structure

- A time-domain pulse-acquire version of the EasySpin biaryl test file,
- with acknowledgements to Stefan Stoll.
- The Spinach simulation is run using explicit time propagation in Liou-
- ville space. Full symmetry treatment using the S2xS2xS2xS2xS2xS2 group
- direct product is performed.
- Calculation time: seconds
- Magnet induction
- Isotope list
- Basis set
- Zeeman interactions
- Spin-spin couplings
- Relaxation theory

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
