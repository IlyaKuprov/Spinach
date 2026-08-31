# examples/nmr_liquids/pansy_triple_ch.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/pansy_triple_ch.m`
- Signature: `pansy_triple_ch()`
- Total lines: 111

## Purpose

Triple-channel PANSY experiment on glycine with natural content of 13C isotope. Coordinates, shieldings, and J- couplings computed with DFT. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Read the spin system properties (vacuum DFT calculation); implemented by `options.min_j=3.0; options.no_xyz=1`.
- Lines 17-18: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 26-27: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 29-30: Sequence parameters; implemented by `parameters.sweep=[1700 4500 8000]`.
- Lines 37-38: Preallocate the answer; implemented by `spec_a=zeros(parameters.zerofill(1),parameters.zerofill(1),'like',1i)`.
- Lines 42-43: Generate isotopomers; implemented by `subsys=dilute(spin_system,'13C')`.
- Lines 45-46: Loop over isotopomers; implemented by `for n=1:numel(subsys)`.
- Lines 48-49: Generate isotopomers; implemented by `subsubsys=dilute(subsys{n,:},'15N')`.
- Lines 51-52: Loop over isotopomers; implemented by `for k=1:numel(subsubsys)`.
- Lines 54-55: Build the basis; implemented by `subsystem=basis(subsubsys{k},bas)`.
- Lines 57-58: Simulation; implemented by `fid=liquid(subsystem,@pansy_triple,parameters,'nmr')`.
- Lines 60-61: Apodisation; implemented by `fid.aa=apodisation(spin_system,fid.aa,{{'sqcos'},{'sqcos'}})`.
- Lines 65-67: Fourier transforms; implemented by `spec_a=spec_a+fftshift(fft2(fid.aa,parameters.zerofill(1), parameters.zerofill(1)))`.
- Lines 77-78: Store the original spectral axes; implemented by `plot_offset=parameters.offset`.
- Lines 83-84: Plotting -1H-15N block; implemented by `kfigure(); scale_figure([3.0 1.5]); subplot(1,3,3)`.
- Lines 92-93: Plotting -1H-13C block; implemented by `subplot(1,3,2)`.
- Lines 101-102: Plotting -1H-1H block; implemented by `subplot(1,3,1)`.

### Control flow inferred from the code

- Line 46: `for` loop over `n=1:numel(subsys)`.
- Line 52: `for` loop over `k=1:numel(subsubsys)`.

### Key state/data transformations

- Lines 13: computes `options.min_j` using `options.min_j=3.0; options.no_xyz=1`.
- Lines 14-15: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/glycine.log'), {{'H','1H'},{'C','13C'},{'N','15N'}},[31.5 189.2 400.3],options)`.
- Lines 18: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 23: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 24: computes `bas.space_level` using `bas.space_level=1`.
- Lines 27: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 30: computes `parameters.sweep` using `parameters.sweep=[1700 4500 8000]`.
- Lines 31: computes `parameters.offset` using `parameters.offset=[900 3200 4000]`.
- Lines 32: computes `parameters.npoints` using `parameters.npoints=[128 128 128]`.
- Lines 33: computes `parameters.zerofill` using `parameters.zerofill=[512 512 512]`.
- Lines 34: computes `parameters.spins` using `parameters.spins={'1H','13C','15N'}`.
- Lines 35: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 38: computes `spec_a` using `spec_a=zeros(parameters.zerofill(1),parameters.zerofill(1),'like',1i)`.
- Lines 39: computes `spec_b` using `spec_b=zeros(parameters.zerofill(2),parameters.zerofill(1),'like',1i)`.
- Lines 40: computes `spec_c` using `spec_c=zeros(parameters.zerofill(3),parameters.zerofill(1),'like',1i)`.
- Lines 43: computes `subsys` using `subsys=dilute(spin_system,'13C')`.

## Implementation structure

- Triple-channel PANSY experiment on glycine with natural
- content of 13C isotope. Coordinates, shieldings, and J-
- couplings computed with DFT.
- Calculation time: seconds
- Read the spin system properties (vacuum DFT calculation)
- Magnet field
- Basis set
- Spinach housekeeping
- Sequence parameters
- Preallocate the answer
- Generate isotopomers
- Loop over isotopomers

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `subplot()`, `plot_offset()`, `plot_sweep()`, `plot_zerofill()`, `plot_spins()`.
