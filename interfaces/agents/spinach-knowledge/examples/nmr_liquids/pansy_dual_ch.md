# examples/nmr_liquids/pansy_dual_ch.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/pansy_dual_ch.m`
- Signature: `pansy_dual_ch()`
- Total lines: 79

## Purpose

PANSY-COSY spectra of camphor with natural content of 13C isotope. Coordinates, shieldings, and J-couplings computed with DFT. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Spin system properties (vacuum DFT calculation); implemented by `options.min_j=3.0; options.no_xyz=1`.
- Lines 16-17: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: Sequence parameters; implemented by `parameters.sweep=[1800 9000]`.
- Lines 36-37: Generate isotopomers; implemented by `subsystems=dilute(spin_system,'13C')`.
- Lines 39-40: Preallocate the answer; implemented by `spec_a=zeros(parameters.zerofill(1),parameters.zerofill(1),'like',1i)`.
- Lines 43-44: Loop over isotopomers; implemented by `parfor n=1:numel(subsystems)`.
- Lines 46-47: Build the basis; implemented by `subsystem=basis(subsystems{n},bas)`.
- Lines 49-50: Simulation; implemented by `fid=liquid(subsystem,@pansy_cosy,parameters,'nmr')`.
- Lines 52-53: Apodisation; implemented by `fid.aa=apodisation(spin_system,fid.aa,{{'sqcos'},{'sqcos'}})`.
- Lines 56-58: Fourier transforms; implemented by `spec_a=spec_a+fftshift(fft2(fid.aa,parameters.zerofill(1), parameters.zerofill(1)))`.
- Lines 64-65: Plotting -heteronuclear side; implemented by `kfigure(); scale_figure([2.0 1.5]); subplot(1,2,2)`.
- Lines 69-70: Plotting -homonuclear side; implemented by `subplot(1,2,1)`.

### Control flow inferred from the code

- Line 44: `parfor` loop over `n=1:numel(subsystems)`.

### Key state/data transformations

- Lines 12: computes `options.min_j` using `options.min_j=3.0; options.no_xyz=1`.
- Lines 13-14: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/camphor.log'), {{'H','1H'},{'C','13C'}},[31.5 189.2],options)`.
- Lines 17: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 22: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 23: computes `bas.space_level` using `bas.space_level=1`.
- Lines 26: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 29: computes `parameters.sweep` using `parameters.sweep=[1800 9000]`.
- Lines 30: computes `parameters.offset` using `parameters.offset=[900 4500]`.
- Lines 31: computes `parameters.npoints` using `parameters.npoints=[128 128]`.
- Lines 32: computes `parameters.zerofill` using `parameters.zerofill=[512 512]`.
- Lines 33: computes `parameters.spins` using `parameters.spins={'1H','13C'}`.
- Lines 34: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 37: computes `subsystems` using `subsystems=dilute(spin_system,'13C')`.
- Lines 40: computes `spec_a` using `spec_a=zeros(parameters.zerofill(1),parameters.zerofill(1),'like',1i)`.
- Lines 41: computes `spec_b` using `spec_b=zeros(parameters.zerofill(2),parameters.zerofill(1),'like',1i)`.
- Lines 47: computes `subsystem` using `subsystem=basis(subsystems{n},bas)`.

## Implementation structure

- PANSY-COSY spectra of camphor with natural content of 13C isotope.
- Coordinates, shieldings, and J-couplings computed with DFT.
- Calculation time: seconds
- Spin system properties (vacuum DFT calculation)
- Magnet field
- Basis set
- Spinach housekeeping
- Sequence parameters
- Generate isotopomers
- Preallocate the answer
- Loop over isotopomers
- Build the basis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `subplot()`, `plot_2d()`.
