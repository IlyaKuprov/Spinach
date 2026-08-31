# examples/nmr_liquids/hmbc_sucrose.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/hmbc_sucrose.m`
- Signature: `hmbc_sucrose()`
- Total lines: 80

## Purpose

HMBC spectrum of sucrose with natural content of 13C isotope (magnetic parameters computed with DFT). Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Spin system properties (vacuum DFT calculation); implemented by `options.min_j=3.0; options.no_xyz=1`.
- Lines 16-17: Set the isotropic parts of shielding tensors to experimental values; implemented by `spin_numbers=[1:19 24:30]`.
- Lines 24-25: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 27-28: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 31-32: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: Sequence parameters; implemented by `parameters.J=140`.
- Lines 47-48: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 50-51: Generate isotopomers; implemented by `subsystems=dilute(spin_system,'13C')`.
- Lines 53-55: Preallocate the answer; implemented by `spectrum=zeros(parameters.zerofill(2), parameters.zerofill(1),'like',1i)`.
- Lines 57-58: Loop over isotopomers; implemented by `parfor n=1:numel(subsystems)`.
- Lines 60-61: Build the basis; implemented by `subsystem=basis(subsystems{n},bas)`.
- Lines 63-64: Simulation; implemented by `fid=liquid(subsystem,@hmbc,parameters,'nmr')`.
- Lines 66-67: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'cos'},{'cos'}})`.
- Lines 69-71: Fourier transform; implemented by `spectrum=spectrum+fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 74-75: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Control flow inferred from the code

- Line 58: `parfor` loop over `n=1:numel(subsystems)`.

### Key state/data transformations

- Lines 12: computes `options.min_j` using `options.min_j=3.0; options.no_xyz=1`.
- Lines 13-14: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'), {{'H','1H'},{'C','13C'}},[31.8 182.1],options)`.
- Lines 17: computes `spin_numbers` using `spin_numbers=[1:19 24:30]`.
- Lines 18-21: computes `new_shifts` using `new_shifts=[ 94.5 73.4 74.9 71.5 74.7 62.4 63.6 106.0 78.7 76.3 83.7 64.7 5.49 3.63 3.83 3.54 3.90 3.90 3.90 3.75 3.75 4.29 4.12 3.96 3.90 3.90]`.
- Lines 22: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,spin_numbers,new_shifts)`.
- Lines 25: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 28: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 29: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 32: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 33: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 34: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 35: computes `bas.space_level` using `bas.space_level=1`.
- Lines 38: computes `parameters.J` using `parameters.J=140`.
- Lines 39: computes `parameters.delta_b` using `parameters.delta_b=60e-3`.
- Lines 40: computes `parameters.sweep` using `parameters.sweep=[6000 2500]`.
- Lines 41: computes `parameters.offset` using `parameters.offset=[5000 900]`.
- Lines 42: computes `parameters.npoints` using `parameters.npoints=[128 128]`.
- Lines 43: computes `parameters.zerofill` using `parameters.zerofill=[512 512]`.

## Implementation structure

- HMBC spectrum of sucrose with natural content of 13C isotope
- (magnetic parameters computed with DFT).
- Calculation time: seconds
- Spin system properties (vacuum DFT calculation)
- Set the isotropic parts of shielding tensors to experimental values
- Magnet field
- Algorithmic options
- Basis set
- Sequence parameters
- Create the spin system structure
- Generate isotopomers
- Preallocate the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
