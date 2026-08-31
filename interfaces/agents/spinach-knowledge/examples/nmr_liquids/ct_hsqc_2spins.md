# examples/nmr_liquids/ct_hsqc_2spins.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/ct_hsqc_2spins.m`
- Signature: `ct_hsqc_2spins()`
- Total lines: 73

## Purpose

CT HSQC spectrum of 2 spin system Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Spin system; implemented by `sys.isotopes={'13C','1H'}`.
- Lines 13-14: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Sequence parameters; implemented by `parameters.J=140`.
- Lines 33-34: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Generate isotopomers; implemented by `subsystems=dilute(spin_system,'13C')`.
- Lines 39-41: Preallocate the answer; implemented by `spectrum=zeros(parameters.zerofill(2), parameters.zerofill(1),'like',1i)`.
- Lines 43-44: Loop over isotopomers; implemented by `parfor n=1:numel(subsystems)`.
- Lines 46-47: Build the basis; implemented by `subsystem=basis(subsystems{n},bas)`.
- Lines 49-50: Simulation; implemented by `fid=liquid(subsystem,@ct_hsqc,parameters,'nmr')`.
- Lines 52-53: Apodisation; implemented by `fid.pos=apodisation(spin_system,fid.pos,{{'sqcos'},{'sqcos'}})`.
- Lines 56-57: F2 Fourier transform; implemented by `f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1)`.
- Lines 60-61: Form States signal; implemented by `fid=f1_pos+conj(f1_neg)`.
- Lines 63-64: F1 Fourier transform; implemented by `spectrum=spectrum+fftshift(fft(fid,parameters.zerofill(1),2),2)`.
- Lines 68-69: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Control flow inferred from the code

- Line 44: `parfor` loop over `n=1:numel(subsystems)`.

### Key state/data transformations

- Lines 11: computes `sys.isotopes` using `sys.isotopes={'13C','1H'}`.
- Lines 14: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={50.00 3.00}`.
- Lines 16: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=140.0`.
- Lines 17: computes `inter.coupling.scalar{2,2}` using `inter.coupling.scalar{2,2}=0`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `parameters.J` using `parameters.J=140`.
- Lines 25: computes `parameters.sweep` using `parameters.sweep=[2500 950]`.
- Lines 26: computes `parameters.offset` using `parameters.offset=[3000 600]`.
- Lines 27: computes `parameters.npoints` using `parameters.npoints=[128 128]`.
- Lines 28: computes `parameters.zerofill` using `parameters.zerofill=[512 512]`.
- Lines 29: computes `parameters.spins` using `parameters.spins={'13C','1H'}`.
- Lines 30: computes `parameters.decouple_f2` using `parameters.decouple_f2={'13C'}`.
- Lines 31: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `subsystems` using `subsystems=dilute(spin_system,'13C')`.
- Lines 40-41: computes `spectrum` using `spectrum=zeros(parameters.zerofill(2), parameters.zerofill(1),'like',1i)`.

## Implementation structure

- CT HSQC spectrum of 2 spin system
- Calculation time: seconds
- Spin system
- Magnet field
- Basis set
- Sequence parameters
- Create the spin system structure
- Generate isotopomers
- Preallocate the answer
- Loop over isotopomers
- Build the basis
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
