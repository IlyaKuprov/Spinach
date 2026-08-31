# examples/nmr_liquids/clip_hsqc_sucrose.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/clip_hsqc_sucrose.m`
- Signature: `clip_hsqc_sucrose()`
- Total lines: 91

## Purpose

CLIP-HSQC spectrum of sucrose with natural content of 13C isotope. Coordinates, shielding anisotropies and J-couplings computed with DFT, isotropic chemical shifts taken from experimental data. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Read the spin system properties (vacuum DFT calculation); implemented by `options.min_j=3.0; options.no_xyz=0`.
- Lines 17-18: Set the isotropic parts of shielding tensors to experimental values; implemented by `spin_numbers=[1:19 24:30]`.
- Lines 25-26: Set the field strength; implemented by `sys.magnet=14.1`.
- Lines 28-29: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 38-39: Sequence parameters; implemented by `parameters.sweep=[8000 2000]`.
- Lines 47-48: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 50-51: Remove fast exchanging and uncoupled spins from the simulation; implemented by `spin_system=kill_spin(spin_system,[20,21,22,23,31,32,33,34])`.
- Lines 53-54: Generate isotopomers; implemented by `subsystems=dilute(spin_system,'13C')`.
- Lines 56-58: Preallocate the result; implemented by `spectrum=zeros(parameters.zerofill(2), parameters.zerofill(1),'like',1i)`.
- Lines 60-61: Run the CLIP-HSQC simulation; implemented by `parfor n=1:numel(subsystems)`.
- Lines 63-64: Build the basis; implemented by `subsystem=basis(subsystems{n},bas)`.
- Lines 66-67: Run simulation; implemented by `fid=liquid(subsystem,@clip_hsqc,parameters,'nmr')`.
- Lines 69-70: Apodisation; implemented by `P_term=apodisation(spin_system,fid.pos,{{'sqcos'},{'sqcos'}})`.
- Lines 73-74: F2 Fourier transform (directly detected dimension); implemented by `f1_P=fftshift(fft(P_term,parameters.zerofill(2),1),1)`.
- Lines 77-78: Form States signal; implemented by `f1_tot=f1_P+conj(f1_N)`.
- Lines 80-81: F1 Fourier transform (indirectly detected dimension); implemented by `spectrum=spectrum+fftshift(fft(f1_tot,parameters.zerofill(1),2),2)`.
- Lines 85-86: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Control flow inferred from the code

- Line 61: `parfor` loop over `n=1:numel(subsystems)`.

### Key state/data transformations

- Lines 13: computes `options.min_j` using `options.min_j=3.0; options.no_xyz=0`.
- Lines 14-15: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'), {{'H','1H'},{'C','13C'}},[31.8 182.1],options)`.
- Lines 18: computes `spin_numbers` using `spin_numbers=[1:19 24:30]`.
- Lines 19-22: computes `new_shifts` using `new_shifts=[ 94.5 73.4 74.9 71.5 74.7 62.4 63.6 106.0 78.7 76.3 83.7 64.7 5.49 3.63 3.83 3.54 3.90 3.90 3.90 3.75 3.75 4.29 4.12 3.96 3.90 3.90]`.
- Lines 23: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,spin_numbers,new_shifts)`.
- Lines 26: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 31: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 32: computes `bas.space_level` using `bas.space_level=1`.
- Lines 35: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 36: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 39: computes `parameters.sweep` using `parameters.sweep=[8000 2000]`.
- Lines 40: computes `parameters.offset` using `parameters.offset=[12000 2700]`.
- Lines 41: computes `parameters.npoints` using `parameters.npoints=[256 256]`.
- Lines 42: computes `parameters.zerofill` using `parameters.zerofill=[512 512]`.
- Lines 43: computes `parameters.spins` using `parameters.spins={'13C','1H'}`.
- Lines 44: computes `parameters.J` using `parameters.J=140`.

## Implementation structure

- CLIP-HSQC spectrum of sucrose with natural content of 13C isotope.
- Coordinates, shielding anisotropies and J-couplings computed with
- DFT, isotropic chemical shifts taken from experimental data.
- Calculation time: minutes
- Read the spin system properties (vacuum DFT calculation)
- Set the isotropic parts of shielding tensors to experimental values
- Set the field strength
- Basis set
- Algorithmic options
- Sequence parameters
- Create the spin system structure
- Remove fast exchanging and uncoupled spins from the simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `kill_spin()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
