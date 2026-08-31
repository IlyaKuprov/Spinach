# examples/nmr_liquids/inad_cyprinol.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/inad_cyprinol.m`
- Signature: `inad_cyprinol()`
- Total lines: 80

## Purpose

INADEQUATE spectrum of cyprinol. The sequence selects double- quantum coherence from coupled 13C pairs and converts it back for detection. A parallel sum over isotopomers that have adjacent 13C spins is used. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Spin system -cyprinol; implemented by `[sys,inter]=cyprinol()`.
- Lines 16-17: Magnet field; implemented by `sys.magnet=11.7`.
- Lines 19-20: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 23-24: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: seq parameters; implemented by `parameters.spins={'13C'}`.
- Lines 41-42: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Generate isotopomers; implemented by `subsystems=dilute(spin_system,'13C',2)`.
- Lines 47-48: Preallocate the answer; implemented by `spectrum=zeros([parameters.zerofill(1) 1],'like',1i)`.
- Lines 50-51: Loop over isotopomers; implemented by `parfor n=1:numel(subsystems)`.
- Lines 53-54: Check the J-coupling between the two carbons; implemented by `c_idx=find(cellfun(@(x)strcmp(x,'13C'),subsystems{n}.comp.isotopes))`.
- Lines 57-58: For couplings stronger than 1 Hz; implemented by `if abs(J)>2*pi*1.0`.
- Lines 60-61: Build the basis; implemented by `subsystem=basis(subsystems{n},bas)`.
- Lines 63-64: Simulation; implemented by `fid=liquid(subsystem,@inadequate,parameters,'nmr')`.
- Lines 66-67: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 69-70: Fourier transform; implemented by `spectrum=spectrum+fftshift(fft(fid,parameters.zerofill))`.
- Lines 76-77: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Control flow inferred from the code

- Line 51: `parfor` loop over `n=1:numel(subsystems)`.
- Line 58: conditional branch on `abs(J)>2*pi*1.0`.

### Key state/data transformations

- Lines 14: computes `[sys,inter]` using `[sys,inter]=cyprinol()`.
- Lines 17: computes `sys.magnet` using `sys.magnet=11.7`.
- Lines 20: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 21: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 24: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 26: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 27: computes `bas.space_level` using `bas.space_level=1`.
- Lines 28: computes `bas.level` using `bas.level=4`.
- Lines 31: computes `parameters.spins` using `parameters.spins={'13C'}`.
- Lines 32: computes `parameters.J` using `parameters.J=50`.
- Lines 33: computes `parameters.decouple` using `parameters.decouple={'1H'}`.
- Lines 34: computes `parameters.offset` using `parameters.offset=5000`.
- Lines 35: computes `parameters.sweep` using `parameters.sweep=10000`.
- Lines 36: computes `parameters.npoints` using `parameters.npoints=4096`.
- Lines 37: computes `parameters.zerofill` using `parameters.zerofill=8192`.
- Lines 38: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 39: computes `parameters.invert_axis` using `parameters.invert_axis=1`.

## Implementation structure

- INADEQUATE spectrum of cyprinol. The sequence selects double-
- quantum coherence from coupled 13C pairs and converts it back
- for detection. A parallel sum over isotopomers that have adjacent
- 13C spins is used.
- Calculation time: minutes
- Spin system -cyprinol
- Magnet field
- Algorithmic options
- Basis set
- seq parameters
- Spinach housekeeping
- Generate isotopomers

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `cyprinol()`, `create()`, `dilute()`, `cellfun()`, `strcmp()`, `get_coupling()`, `c_idx()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
