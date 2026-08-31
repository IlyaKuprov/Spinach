# examples/nmr_proteins/noesyhsqc_ubiquitin_prot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/noesyhsqc_ubiquitin_prot.m`
- Signature: `noesyhsqc_ubiquitin_prot()`
- Total lines: 94

## Purpose

1H-1H-15N NOESY-HSQC spectrum of 15N-labelled ubiquitin at 900 MHz with 65 ms mixing time. It is assumed that the protein is not 13C-labelled. Calculation time: hours.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Protein data import; implemented by `options.pdb_mol=1`.
- Lines 18-19: Magnet field; implemented by `sys.magnet=21.1356`.
- Lines 21-22: Tolerances; implemented by `sys.tols.inter_cutoff=2.0`.
- Lines 25-26: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 31-32: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: Algorithmic options; implemented by `sys.enable={'prop_cache','greedy'}`.
- Lines 41-42: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Kill carbons; implemented by `spin_system=kill_spin(spin_system,strcmp('13C',spin_system.comp.isotopes))`.
- Lines 47-48: Build the basis; implemented by `spin_system=basis(spin_system,bas)`.
- Lines 50-51: Sequence parameters; implemented by `parameters.tmix=0.065`.
- Lines 60-61: Simulation; implemented by `fid=liquid(spin_system,@noesyhsqc,parameters,'nmr')`.
- Lines 63-64: Apodisation; implemented by `fid.pos_pos=apodisation(spin_system,fid.pos_pos,{{'sqcos'},{'sqcos'},{'sqcos'}})`.
- Lines 69-70: F3 Fourier transform; implemented by `f3_pos_pos=fftshift(fft(fid.pos_pos,parameters.zerofill(3),3),3)`.
- Lines 75-76: Absorption part of F3 signal; implemented by `f3_pos=f3_pos_pos+conj(f3_neg_neg)`.
- Lines 79-80: F2 Fourier transform; implemented by `f3f2_pos=fftshift(fft(f3_pos,parameters.zerofill(2),2),2)`.
- Lines 83-84: Absorption part of F2 signal; implemented by `f3f2=f3f2_pos+conj(f3f2_neg)`.
- Lines 86-87: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f3f2,parameters.zerofill(1),1),1)`.
- Lines 89-91: Plotting; implemented by `kfigure(); plot_3d(spin_system,-real(spectrum),parameters, 10,[0.01 0.1 0.01 0.1],2,'positive')`.

### Key state/data transformations

- Lines 13: computes `options.pdb_mol` using `options.pdb_mol=1`.
- Lines 14: computes `options.select` using `options.select='all'`.
- Lines 15: computes `options.noshift` using `options.noshift='delete'`.
- Lines 16: computes `[sys,inter]` using `[sys,inter]=protein('1D3Z.pdb','1D3Z.bmrb',options)`.
- Lines 19: computes `sys.magnet` using `sys.magnet=21.1356`.
- Lines 22: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=2.0`.
- Lines 23: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 26: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 27: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 28: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 29: computes `inter.tau_c` using `inter.tau_c={5e-9}`.
- Lines 32: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 33: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 34: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 35: computes `bas.level` using `bas.level=4; bas.space_level=3`.
- Lines 38: computes `sys.enable` using `sys.enable={'prop_cache','greedy'}`.
- Lines 39: computes `sys.disable` using `sys.disable={'asyredf'}`.
- Lines 42: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- 1H-1H-15N NOESY-HSQC spectrum of 15N-labelled ubiquitin at 900
- MHz with 65 ms mixing time. It is assumed that the protein is
- not 13C-labelled.
- Calculation time: hours.
- Protein data import
- Magnet field
- Tolerances
- Relaxation theory
- Basis set
- Algorithmic options
- Create the spin system structure
- Kill carbons

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `protein()`, `create()`, `kill_spin()`, `strcmp()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `plot_3d()`.
