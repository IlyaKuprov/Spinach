# examples/nmr_liquids/inept_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/inept_strychnine.m`
- Signature: `inept_strychnine()`
- Total lines: 66

## Purpose

INEPT experiment on strychnine. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Read the spin system properties; implemented by `[sys,inter]=strychnine({'1H','13C'})`.
- Lines 13-14: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 16-17: Temperature; implemented by `inter.temperature=298`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Sequence parameters; implemented by `parameters.sweep=10000`.
- Lines 33-34: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Generate isotopomers; implemented by `subsystems=dilute(spin_system,'13C')`.
- Lines 39-40: Preallocate the answer; implemented by `spectrum=zeros([parameters.zerofill 1],'like',1i)`.
- Lines 42-43: Loop over isotopomers; implemented by `parfor n=1:numel(subsystems)`.
- Lines 45-46: Build the basis; implemented by `subsystem=basis(subsystems{n},bas)`.
- Lines 48-49: Simulation; implemented by `fid=liquid(subsystem,@inept,parameters,'nmr')`.
- Lines 51-52: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 54-55: Fourier transform; implemented by `spectrum=spectrum+fftshift(fft(fid,parameters.zerofill))`.
- Lines 59-60: Plotting; implemented by `kfigure(); parameters.spins={'13C'}`.

### Control flow inferred from the code

- Line 43: `parfor` loop over `n=1:numel(subsystems)`.

### Key state/data transformations

- Lines 11: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H','13C'})`.
- Lines 14: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 17: computes `inter.temperature` using `inter.temperature=298`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 22: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 23: computes `bas.space_level` using `bas.space_level=1`.
- Lines 26: computes `parameters.sweep` using `parameters.sweep=10000`.
- Lines 27: computes `parameters.offset` using `parameters.offset=[5000 0]`.
- Lines 28: computes `parameters.npoints` using `parameters.npoints=2048`.
- Lines 29: computes `parameters.zerofill` using `parameters.zerofill=8196`.
- Lines 30: computes `parameters.spins` using `parameters.spins={'13C','1H'}`.
- Lines 31: computes `parameters.J` using `parameters.J=150`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `subsystems` using `subsystems=dilute(spin_system,'13C')`.
- Lines 40: computes `spectrum` using `spectrum=zeros([parameters.zerofill 1],'like',1i)`.
- Lines 46: computes `subsystem` using `subsystem=basis(subsystems{n},bas)`.
- Lines 49: computes `fid` using `fid=liquid(subsystem,@inept,parameters,'nmr')`.

## Implementation structure

- INEPT experiment on strychnine.
- Calculation time: minutes
- Read the spin system properties
- Magnet field
- Temperature
- Basis set
- Sequence parameters
- Create the spin system structure
- Generate isotopomers
- Preallocate the answer
- Loop over isotopomers
- Build the basis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
