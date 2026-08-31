# examples/nmr_zerofield/field_drop_acetonitrile.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_zerofield/field_drop_acetonitrile.m`
- Signature: `field_drop_acetonitrile()`
- Total lines: 76

## Purpose

Zero-field NMR spectroscopy -acetonitrile. The simulation proceeds by computing the exact thermal equilibrium state and them propagating it through a time-dependent field drop. Set to reproduce Figure 7 from Calculation time: seconds

## Physical / mathematical content

- Zero- and ultralow-field NMR examples. The main physics is the crossover from Zeeman-dominated spectra to J-dominated spectra, with coherent evolution in near-zero field and detection of low-frequency transitions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnetic field (polariser); implemented by `sys.magnet=2.0`.
- Lines 17-18: Spin system; implemented by `sys.isotopes={'1H','1H','1H','13C','13C','15N'}`.
- Lines 20-21: Interactions; implemented by `inter.coupling.scalar{1,4}=136.200`.
- Lines 35-36: Temperature; implemented by `inter.temperature=298`.
- Lines 38-39: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 44-45: Sequence parameters; implemented by `parameters.sweep=700`.
- Lines 59-60: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 63-64: Simulation; implemented by `fid=liquid(spin_system,@zulf_abrupt,parameters,'labframe')`.
- Lines 66-67: Apodisation; implemented by `fid=apodisation(spin_system,fid-mean(fid),{{'exp',6}})`.
- Lines 69-70: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 72-73: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=2.0`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','13C','13C','15N'}`.
- Lines 21: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=136.200`.
- Lines 22: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=136.200`.
- Lines 23: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=136.200`.
- Lines 24: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}=-9.924`.
- Lines 25: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}=-9.924`.
- Lines 26: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=-9.924`.
- Lines 27: computes `inter.coupling.scalar{1,6}` using `inter.coupling.scalar{1,6}=-1.688`.
- Lines 28: computes `inter.coupling.scalar{2,6}` using `inter.coupling.scalar{2,6}=-1.688`.
- Lines 29: computes `inter.coupling.scalar{3,6}` using `inter.coupling.scalar{3,6}=-1.688`.
- Lines 30: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=57.010`.
- Lines 31: computes `inter.coupling.scalar{4,6}` using `inter.coupling.scalar{4,6}=2.822`.
- Lines 32: computes `inter.coupling.scalar{5,6}` using `inter.coupling.scalar{5,6}=-17.419`.
- Lines 33: computes `inter.coupling.scalar{6,6}` using `inter.coupling.scalar{6,6}=0`.
- Lines 36: computes `inter.temperature` using `inter.temperature=298`.
- Lines 39: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 40: computes `bas.approximation` using `bas.approximation='none'`.

## Implementation structure

- Zero-field NMR spectroscopy -acetonitrile. The simulation
- proceeds by computing the exact thermal equilibrium state
- and them propagating it through a time-dependent field drop.
- Set to reproduce Figure 7 from
- Calculation time: seconds
- Magnetic field (polariser)
- Spin system
- Interactions
- Temperature
- Basis set
- Sequence parameters
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
