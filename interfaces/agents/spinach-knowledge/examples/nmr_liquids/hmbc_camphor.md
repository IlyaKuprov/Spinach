# examples/nmr_liquids/hmbc_camphor.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/hmbc_camphor.m`
- Signature: `hmbc_camphor()`
- Total lines: 74

## Purpose

HMBC spectrum of camphor with natural content of 13C isotope. Coordinates, shielding anisotropies and J-couplings computed witt DFT. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Spin system properties (vacuum DFT calculation); implemented by `options.min_j=3.0; options.no_xyz=0`.
- Lines 17-18: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 20-21: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 24-25: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Sequence parameters; implemented by `parameters.J=140`.
- Lines 40-41: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 43-44: Generate isotopomers; implemented by `subsystems=dilute(spin_system,'13C')`.
- Lines 46-48: Preallocate the answer; implemented by `spectrum=zeros(parameters.zerofill(2), parameters.zerofill(1),'like',1i)`.
- Lines 50-51: Loop over isotopomers; implemented by `parfor n=1:numel(subsystems)`.
- Lines 53-54: Build the basis; implemented by `subsystem=basis(subsystems{n},bas)`.
- Lines 56-57: Simulation; implemented by `fid=liquid(subsystem,@hmbc,parameters,'nmr')`.
- Lines 59-60: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'cos'},{'cos'}})`.
- Lines 62-64: Fourier transform; implemented by `spectrum=spectrum+fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 68-69: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Control flow inferred from the code

- Line 51: `parfor` loop over `n=1:numel(subsystems)`.

### Key state/data transformations

- Lines 13: computes `options.min_j` using `options.min_j=3.0; options.no_xyz=0`.
- Lines 14-15: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/camphor.log'), {{'H','1H'},{'C','13C'}},[31.8-0.35 182.1+7.14],options)`.
- Lines 18: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 21: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 22: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 25: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 27: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 28: computes `bas.space_level` using `bas.space_level=1`.
- Lines 31: computes `parameters.J` using `parameters.J=140`.
- Lines 32: computes `parameters.delta_b` using `parameters.delta_b=60e-3`.
- Lines 33: computes `parameters.sweep` using `parameters.sweep=[40000 1500]`.
- Lines 34: computes `parameters.offset` using `parameters.offset=[18000 900]`.
- Lines 35: computes `parameters.npoints` using `parameters.npoints=[128 128]`.
- Lines 36: computes `parameters.zerofill` using `parameters.zerofill=[512 512]`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'13C','1H'}`.
- Lines 38: computes `parameters.axis_units` using `parameters.axis_units='Hz'`.
- Lines 41: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- HMBC spectrum of camphor with natural content of 13C isotope.
- Coordinates, shielding anisotropies and J-couplings computed
- witt DFT.
- Calculation time: seconds
- Spin system properties (vacuum DFT calculation)
- Magnet field
- Algorithmic options
- Basis set
- Sequence parameters
- Create the spin system structure
- Generate isotopomers
- Preallocate the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
