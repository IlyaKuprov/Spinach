# examples/nmr_proteins/hcch_cosy_gb1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/hcch_cosy_gb1.m`
- Signature: `hcch_cosy_gb1()`
- Total lines: 87

## Purpose

3D HCCH COSY experiment on GB1 protein. Calculation time: hours.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Protein data import; implemented by `options.pdb_mol=1`.
- Lines 16-17: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 19-20: Tolerances; implemented by `sys.tols.inter_cutoff=20.0`.
- Lines 23-24: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 29-30: Algorithmic options; implemented by `sys.enable={'greedy','prop_cache'}`.
- Lines 32-33: Sequence parameters; implemented by `parameters.J_ch=140`.
- Lines 44-45: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Kill nitrogens (not relevant); implemented by `spin_system=kill_spin(spin_system,strcmp('15N',spin_system.comp.isotopes))`.
- Lines 50-51: Build the basis; implemented by `spin_system=basis(spin_system,bas)`.
- Lines 53-54: Simulation; implemented by `fid=liquid(spin_system,@hcch_cosy,parameters,'nmr')`.
- Lines 56-57: Apodisation; implemented by `fid.pos_pos=apodisation(spin_system,fid.pos_pos,{{'sqcos'},{'sqcos'},{'sqcos'}})`.
- Lines 62-63: F3 Fourier transform; implemented by `f3_pos_pos=fftshift(fft(fid.pos_pos,parameters.zerofill(3),3),3)`.
- Lines 68-69: Absorption part of F3 signal; implemented by `f3_pos=f3_pos_pos+conj(f3_neg_neg)`.
- Lines 72-73: F2 Fourier transform; implemented by `f3f2_pos=fftshift(fft(f3_pos,parameters.zerofill(2),2),2)`.
- Lines 76-77: Absorption part of F2 signal; implemented by `f3f2=f3f2_pos+conj(f3f2_neg)`.
- Lines 79-80: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f3f2,parameters.zerofill(1),1),1)`.
- Lines 82-84: Plotting; implemented by `kfigure(); plot_3d(spin_system,imag(spectrum),parameters, 10,[0.05 0.25 0.05 0.25],2,'positive')`.

### Key state/data transformations

- Lines 11: computes `options.pdb_mol` using `options.pdb_mol=1`.
- Lines 12: computes `options.noshift` using `options.noshift='delete'`.
- Lines 13: computes `options.select` using `options.select='all'`.
- Lines 14: computes `[sys,inter]` using `[sys,inter]=protein('2N9K.pdb','2N9K.bmrb',options)`.
- Lines 17: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 20: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=20.0`.
- Lines 21: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 24: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 26: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 27: computes `bas.level` using `bas.level=4; bas.space_level=1`.
- Lines 30: computes `sys.enable` using `sys.enable={'greedy','prop_cache'}`.
- Lines 33: computes `parameters.J_ch` using `parameters.J_ch=140`.
- Lines 34: computes `parameters.J_cc` using `parameters.J_cc=35`.
- Lines 35: computes `parameters.delta` using `parameters.delta=1.1e-3`.
- Lines 36: computes `parameters.sweep` using `parameters.sweep=[6000 13000 6000]`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'1H','13C','1H'}`.
- Lines 38: computes `parameters.offset` using `parameters.offset=[2500 7000 2500]`.

## Implementation structure

- 3D HCCH COSY experiment on GB1 protein.
- Calculation time: hours.
- Protein data import
- Magnet field
- Tolerances
- Basis set
- Algorithmic options
- Sequence parameters
- Create the spin system structure
- Kill nitrogens (not relevant)
- Build the basis
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `protein()`, `create()`, `kill_spin()`, `strcmp()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `plot_3d()`.
