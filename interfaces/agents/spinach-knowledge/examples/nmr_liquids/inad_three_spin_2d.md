# examples/nmr_liquids/inad_three_spin_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/inad_three_spin_2d.m`
- Signature: `inad_three_spin_2d()`
- Total lines: 81

## Purpose

Example of a 2D-INADEQUATE spectrum of a generic three-spin system. Calculation time: seconds. Theresa Hune Christian Griesinger

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnetic field (700 MHz); implemented by `sys.magnet=16.44`.
- Lines 14-15: Generic three-spin system; implemented by `sys.isotopes={'13C','13C','13C'}`.
- Lines 22-23: Formalism and basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 26-27: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Generate isotopomers; implemented by `subsystems=dilute(spin_system,'13C',2)`.
- Lines 33-34: Sequence parameters; implemented by `parameters.offset=17604.78`.
- Lines 43-45: Preallocate the answer; implemented by `spectrum=zeros(parameters.zerofill(2), parameters.zerofill(1))`.
- Lines 47-48: Loop over isotopomers; implemented by `parfor n=1:numel(subsystems)`.
- Lines 50-51: Build the basis; implemented by `subsystem=basis(subsystems{n},bas)`.
- Lines 53-54: Simulation; implemented by `fid=liquid(subsystem,@inadequate_2d,parameters,'nmr')`.
- Lines 56-57: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'cos'},{'cos'}})`.
- Lines 60-61: F2 Fourier transform; implemented by `spec_cos=fftshift(fft(fid.cos,parameters.zerofill(2),1),1)`.
- Lines 64-65: Form States signal; implemented by `spec_states=real(spec_cos)+1i*real(spec_sin)`.
- Lines 67-68: F1 Fourier transform; implemented by `spectrum=spectrum+fftshift(fft(spec_states,parameters.zerofill(1),2),2)`.
- Lines 72-73: Do the plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Control flow inferred from the code

- Line 48: `parfor` loop over `n=1:numel(subsystems)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=16.44`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'13C','13C','13C'}`.
- Lines 16: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={10 30 70}`.
- Lines 17: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=20`.
- Lines 18: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=60`.
- Lines 19: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=0`.
- Lines 20: computes `inter.coupling.scalar{3,3}` using `inter.coupling.scalar{3,3}=0`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `subsystems` using `subsystems=dilute(spin_system,'13C',2)`.
- Lines 34: computes `parameters.offset` using `parameters.offset=17604.78`.
- Lines 35: computes `parameters.J` using `parameters.J=50`.
- Lines 36: computes `parameters.spins` using `parameters.spins={'13C'}`.
- Lines 37: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 38: computes `parameters.sweep` using `parameters.sweep=[35213.086 34722.223]`.
- Lines 39: computes `parameters.npoints` using `parameters.npoints=[128 2048]`.
- Lines 40: computes `parameters.zerofill` using `parameters.zerofill=[512 8192]`.

## Implementation structure

- Example of a 2D-INADEQUATE spectrum of a
- generic three-spin system.
- Calculation time: seconds.
- Theresa Hune
- Christian Griesinger
- Magnetic field (700 MHz)
- Generic three-spin system
- Formalism and basis set
- Spinach housekeeping
- Generate isotopomers
- Sequence parameters
- Preallocate the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `dilute()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`, `kylabel()`.
