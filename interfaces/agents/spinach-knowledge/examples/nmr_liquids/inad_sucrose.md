# examples/nmr_liquids/inad_sucrose.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/inad_sucrose.m`
- Signature: `inad_sucrose()`
- Total lines: 89

## Purpose

INADEQUATE spectrum of sucrose. The sequence selects double- quantum coherence from coupled 13C pairs and converts it back for detection. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Spin system properties (vacuum DFT calculation); implemented by `options.min_j=5.0; options.no_xyz=1`.
- Lines 17-18: Set the isotropic parts of shielding tensors to experimental values; implemented by `spin_numbers=[1:19 24:30]`.
- Lines 25-26: Magnet field; implemented by `sys.magnet=11.7`.
- Lines 28-29: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 39-40: Sequence parameters; implemented by `parameters.spins={'13C'}`.
- Lines 50-51: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 53-54: Generate isotopomers; implemented by `subsystems=dilute(spin_system,'13C',2)`.
- Lines 56-57: Preallocate the answer; implemented by `spectrum=zeros([parameters.zerofill(1) 1],'like',1i)`.
- Lines 59-60: Loop over isotopomers; implemented by `parfor n=1:numel(subsystems)`.
- Lines 62-63: Check the J-coupling between the two carbons; implemented by `c_idx=find(cellfun(@(x)strcmp(x,'13C'),subsystems{n}.comp.isotopes))`.
- Lines 66-67: For couplings stronger than 1 Hz; implemented by `if abs(J)>2*pi*1.0`.
- Lines 69-70: Build the basis; implemented by `subsystem=basis(subsystems{n},bas)`.
- Lines 72-73: Simulation; implemented by `fid=liquid(subsystem,@inadequate,parameters,'nmr')`.
- Lines 75-76: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 78-79: Fourier transform; implemented by `spectrum=spectrum+fftshift(fft(fid,parameters.zerofill))`.
- Lines 85-86: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Control flow inferred from the code

- Line 60: `parfor` loop over `n=1:numel(subsystems)`.
- Line 67: conditional branch on `abs(J)>2*pi*1.0`.

### Key state/data transformations

- Lines 13: computes `options.min_j` using `options.min_j=5.0; options.no_xyz=1`.
- Lines 14-15: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'), {{'H','1H'},{'C','13C'}},[31.8 189.7],options)`.
- Lines 18: computes `spin_numbers` using `spin_numbers=[1:19 24:30]`.
- Lines 19-22: computes `new_shifts` using `new_shifts=[ 94.5 73.4 74.9 71.5 74.7 62.4 63.6 106.0 78.7 76.3 83.7 64.7 5.49 3.63 3.83 3.54 3.90 3.90 3.90 3.75 3.75 4.29 4.12 3.96 3.90 3.90]`.
- Lines 23: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,spin_numbers,new_shifts)`.
- Lines 26: computes `sys.magnet` using `sys.magnet=11.7`.
- Lines 29: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 30: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 35: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 36: computes `bas.space_level` using `bas.space_level=1`.
- Lines 37: computes `bas.level` using `bas.level=4`.
- Lines 40: computes `parameters.spins` using `parameters.spins={'13C'}`.
- Lines 41: computes `parameters.J` using `parameters.J=50`.
- Lines 42: computes `parameters.decouple` using `parameters.decouple={'1H'}`.
- Lines 43: computes `parameters.offset` using `parameters.offset=10000`.
- Lines 44: computes `parameters.sweep` using `parameters.sweep=8000`.

## Implementation structure

- INADEQUATE spectrum of sucrose. The sequence selects double-
- quantum coherence from coupled 13C pairs and converts it back
- for detection.
- Calculation time: minutes
- Spin system properties (vacuum DFT calculation)
- Set the isotropic parts of shielding tensors to experimental values
- Magnet field
- Algorithmic options
- Basis set
- Sequence parameters
- Spinach housekeeping
- Generate isotopomers

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `dilute()`, `cellfun()`, `strcmp()`, `get_coupling()`, `c_idx()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
