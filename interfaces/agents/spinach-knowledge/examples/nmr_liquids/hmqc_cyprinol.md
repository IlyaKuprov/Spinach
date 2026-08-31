# examples/nmr_liquids/hmqc_cyprinol.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/hmqc_cyprinol.m`
- Signature: `hmqc_cyprinol()`
- Total lines: 73

## Purpose

HMQC spectrum of cyprinol with natural abundance of 13C isotope. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Read the spin system properties; implemented by `[sys,inter,bas]=cyprinol()`.
- Lines 13-14: Magnet field; implemented by `sys.magnet=11.7`.
- Lines 16-17: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 28-29: Sequence parameters; implemented by `parameters.J=150`.
- Lines 39-40: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Generate isotopomers; implemented by `subsystems=dilute(spin_system,'13C')`.
- Lines 45-47: Preallocate the answer; implemented by `spectrum=zeros(parameters.zerofill(2), parameters.zerofill(1),'like',1i)`.
- Lines 49-50: Loop over isotopomers; implemented by `parfor n=1:numel(subsystems)`.
- Lines 52-53: Build the basis; implemented by `subsystem=basis(subsystems{n},bas)`.
- Lines 55-56: Simulation; implemented by `fid=liquid(subsystem,@hmqc,parameters,'nmr')`.
- Lines 58-59: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'cos'},{'cos'}})`.
- Lines 61-63: Fourier transform; implemented by `spectrum=spectrum+fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 67-68: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Control flow inferred from the code

- Line 50: `parfor` loop over `n=1:numel(subsystems)`.

### Key state/data transformations

- Lines 11: computes `[sys,inter,bas]` using `[sys,inter,bas]=cyprinol()`.
- Lines 14: computes `sys.magnet` using `sys.magnet=11.7`.
- Lines 17: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 18: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 19: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=5.0`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 24: computes `bas.level` using `bas.level=3`.
- Lines 25: computes `bas.space_level` using `bas.space_level=1`.
- Lines 26: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 29: computes `parameters.J` using `parameters.J=150`.
- Lines 30: computes `parameters.sweep` using `parameters.sweep=[12000 2500]`.
- Lines 31: computes `parameters.offset` using `parameters.offset=[5000 1250]`.
- Lines 32: computes `parameters.npoints` using `parameters.npoints=[128 128]`.
- Lines 33: computes `parameters.zerofill` using `parameters.zerofill=[512 512]`.
- Lines 34: computes `parameters.spins` using `parameters.spins={'13C','1H'}`.
- Lines 35: computes `parameters.decouple_f1` using `parameters.decouple_f1={'1H'}`.
- Lines 36: computes `parameters.decouple_f2` using `parameters.decouple_f2={'13C'}`.

## Implementation structure

- HMQC spectrum of cyprinol with natural abundance of 13C isotope.
- Calculation time: seconds
- Read the spin system properties
- Magnet field
- Algorithmic options
- Basis set
- Sequence parameters
- Create the spin system structure
- Generate isotopomers
- Preallocate the answer
- Loop over isotopomers
- Build the basis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cyprinol()`, `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
