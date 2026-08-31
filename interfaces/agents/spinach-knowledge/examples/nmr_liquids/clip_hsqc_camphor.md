# examples/nmr_liquids/clip_hsqc_camphor.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/clip_hsqc_camphor.m`
- Signature: `clip_hsqc_camphor()`
- Total lines: 87

## Purpose

CLIP-HSQC spectrum of camphor with natural content of 13C isotope. Coordinates, shielding anisotropies and J-couplings computed with DFT, isotropic chemical shifts taken from experimental data. Calculation time: minutes.

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Spin system properties (vacuum DFT calculation); implemented by `options.min_j=3.0; options.no_xyz=0`.
- Lines 16-17: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 19-25: Isotropic shift components come from the experiment; implemented by `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1:26, [ 29.70 26.80 44.20 59.70 42.80 218.1 49.70 17.80 17.30 7.80 1.35 1.67 1.97 1.31 1.96 0.85 0.85 0.85 0.99 0.99 0.…`.
- Lines 27-28: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 33-34: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 37-38: Sequence parameters; implemented by `parameters.J=140`.
- Lines 46-47: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 49-50: Generate isotopomers; implemented by `subsystems=dilute(spin_system,'13C')`.
- Lines 52-54: Preallocate the answer; implemented by `spectrum=zeros(parameters.zerofill(2), parameters.zerofill(1),'like',1i)`.
- Lines 56-57: Loop over isotopomers; implemented by `parfor n=1:numel(subsystems)`.
- Lines 59-60: Build the basis; implemented by `subsystem=basis(subsystems{n},bas)`.
- Lines 62-63: Simulation; implemented by `fid=liquid(subsystem,@clip_hsqc,parameters,'nmr')`.
- Lines 65-66: Apodisation; implemented by `fid.pos=apodisation(spin_system,fid.pos,{{'sqcos'},{'sqcos'}})`.
- Lines 69-70: F2 Fourier transform; implemented by `f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1)`.
- Lines 73-74: Form States signal; implemented by `fid=f1_pos+conj(f1_neg)`.
- Lines 76-77: F1 Fourier transform; implemented by `spectrum=spectrum+fftshift(fft(fid,parameters.zerofill(1),2),2)`.
- Lines 81-82: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Control flow inferred from the code

- Line 57: `parfor` loop over `n=1:numel(subsystems)`.

### Key state/data transformations

- Lines 13: computes `options.min_j` using `options.min_j=3.0; options.no_xyz=0`.
- Lines 14-15: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/camphor.log'), {{'H','1H'},{'C','13C'}},[31.8 182.1],options)`.
- Lines 17: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 20-25: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1:26, [ 29.70 26.80 44.20 59.70 42.80 218.1 49.70 17.80 17.30 7.80 1.35 1.67 1.97 1.31 1.96 0.85 0.85 0.85 0.99 0.99 0.…`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 30: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 31: computes `bas.space_level` using `bas.space_level=1`.
- Lines 34: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 35: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 38: computes `parameters.J` using `parameters.J=140`.
- Lines 39: computes `parameters.sweep` using `parameters.sweep=[8000 1500]`.
- Lines 40: computes `parameters.offset` using `parameters.offset=[4000 1000]`.
- Lines 41: computes `parameters.npoints` using `parameters.npoints=[128 128]`.
- Lines 42: computes `parameters.zerofill` using `parameters.zerofill=[512 512]`.
- Lines 43: computes `parameters.spins` using `parameters.spins={'13C','1H'}`.
- Lines 44: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 47: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- CLIP-HSQC spectrum of camphor with natural content of 13C isotope.
- Coordinates, shielding anisotropies and J-couplings computed with
- DFT, isotropic chemical shifts taken from experimental data.
- Calculation time: minutes.
- Spin system properties (vacuum DFT calculation)
- Magnet field
- Isotropic shift components come from the experiment
- Basis set
- Algorithmic options
- Sequence parameters
- Create the spin system structure
- Generate isotopomers

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
