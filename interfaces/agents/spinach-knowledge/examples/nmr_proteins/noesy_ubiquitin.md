# examples/nmr_proteins/noesy_ubiquitin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/noesy_ubiquitin.m`
- Signature: `noesy_ubiquitin()`
- Total lines: 82

## Purpose

1H-1H NOESY spectrum of ubiquitin with 65 ms mixing time. It is assumed that the protein is not 13C-or 15N-labelled. Calculation time: hours.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Protein data import; implemented by `options.pdb_mol=1`.
- Lines 17-18: Magnet field; implemented by `sys.magnet=21.1356`.
- Lines 20-21: Tolerances; implemented by `sys.tols.inter_cutoff=2.0`.
- Lines 24-25: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Algorithmic options; implemented by `sys.enable={'prop_cache','greedy'}`.
- Lines 39-40: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Kill carbons and nitrogens (protein assumed unlabelled); implemented by `spin_system=kill_spin(spin_system,strcmp('13C',spin_system.comp.isotopes))`.
- Lines 46-47: Build the basis; implemented by `spin_system=basis(spin_system,bas)`.
- Lines 49-50: Sequence parameters; implemented by `parameters.tmix=0.065`.
- Lines 59-60: Simulation; implemented by `fid=liquid(spin_system,@noesy,parameters,'nmr')`.
- Lines 62-63: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 66-67: F2 Fourier transform; implemented by `f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.
- Lines 70-71: States signal; implemented by `f1_states=f1_cos-1i*f1_sin`.
- Lines 73-74: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 76-77: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 12: computes `options.pdb_mol` using `options.pdb_mol=1`.
- Lines 13: computes `options.select` using `options.select='all'`.
- Lines 14: computes `options.noshift` using `options.noshift='delete'`.
- Lines 15: computes `[sys,inter]` using `[sys,inter]=protein('1D3Z.pdb','1D3Z.bmrb',options)`.
- Lines 18: computes `sys.magnet` using `sys.magnet=21.1356`.
- Lines 21: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=2.0`.
- Lines 22: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 27: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 28: computes `inter.tau_c` using `inter.tau_c={5e-9}`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 33: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 34: computes `bas.level` using `bas.level=4; bas.space_level=3`.
- Lines 37: computes `sys.enable` using `sys.enable={'prop_cache','greedy'}`.
- Lines 40: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 50: computes `parameters.tmix` using `parameters.tmix=0.065`.

## Implementation structure

- 1H-1H NOESY spectrum of ubiquitin with 65 ms mixing time. It is
- assumed that the protein is not 13C-or 15N-labelled.
- Calculation time: hours.
- Protein data import
- Magnet field
- Tolerances
- Relaxation theory
- Basis set
- Algorithmic options
- Create the spin system structure
- Kill carbons and nitrogens (protein assumed unlabelled)
- Build the basis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `protein()`, `create()`, `kill_spin()`, `strcmp()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
