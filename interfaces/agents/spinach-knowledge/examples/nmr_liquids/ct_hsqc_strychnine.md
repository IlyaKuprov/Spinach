# examples/nmr_liquids/ct_hsqc_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/ct_hsqc_strychnine.m`
- Signature: `ct_hsqc_strychnine()`
- Total lines: 76

## Purpose

CT HSQC spectrum of strychnine with natural content of 13C isotope. Calculation time: hours.

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Read the spin system properties; implemented by `[sys,inter]=strychnine({'13C','1H'})`.
- Lines 13-14: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 16-17: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 26-27: Sequence parameters; implemented by `parameters.J=140`.
- Lines 36-37: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Generate isotopomers; implemented by `subsystems=dilute(spin_system,'13C')`.
- Lines 42-44: Preallocate the answer; implemented by `spectrum=zeros(parameters.zerofill(2), parameters.zerofill(1),'like',1i)`.
- Lines 46-47: Loop over isotopomers; implemented by `parfor n=1:numel(subsystems)`.
- Lines 49-50: Build the basis; implemented by `subsystem=basis(subsystems{n},bas)`.
- Lines 52-53: Simulation; implemented by `fid=liquid(subsystem,@ct_hsqc,parameters,'nmr')`.
- Lines 55-56: Apodisation; implemented by `fid.pos=apodisation(spin_system,fid.pos,{{'sqcos'},{'sqcos'}})`.
- Lines 59-60: F2 Fourier transform; implemented by `f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1)`.
- Lines 63-64: Form States signal; implemented by `fid=f1_pos+conj(f1_neg)`.
- Lines 66-67: F1 Fourier transform; implemented by `spectrum=spectrum+fftshift(fft(fid,parameters.zerofill(1),2),2)`.
- Lines 71-72: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Control flow inferred from the code

- Line 47: `parfor` loop over `n=1:numel(subsystems)`.

### Key state/data transformations

- Lines 11: computes `[sys,inter]` using `[sys,inter]=strychnine({'13C','1H'})`.
- Lines 14: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 17: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 18: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 23: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 24: computes `bas.space_level` using `bas.space_level=1`.
- Lines 27: computes `parameters.J` using `parameters.J=140`.
- Lines 28: computes `parameters.sweep` using `parameters.sweep=[10000 3000]`.
- Lines 29: computes `parameters.offset` using `parameters.offset=[4000 1000]`.
- Lines 30: computes `parameters.npoints` using `parameters.npoints=[256 256]`.
- Lines 31: computes `parameters.zerofill` using `parameters.zerofill=[512 512]`.
- Lines 32: computes `parameters.spins` using `parameters.spins={'13C','1H'}`.
- Lines 33: computes `parameters.decouple_f2` using `parameters.decouple_f2={'13C'}`.
- Lines 34: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `subsystems` using `subsystems=dilute(spin_system,'13C')`.

## Implementation structure

- CT HSQC spectrum of strychnine with natural content of 13C isotope.
- Calculation time: hours.
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

- Called routines detected from the main body: `strychnine()`, `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
