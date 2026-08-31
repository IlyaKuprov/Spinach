# examples/nmr_proteins/methyl_noesy_gb1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/methyl_noesy_gb1.m`
- Signature: `methyl_noesy_gb1()`
- Total lines: 84

## Purpose

1H-1H NOESY spectrum of GB1 with everything deuterated except methyl groups. Deuteria are kept in the spin system because they are a part of the coupling network; methyl group rotati- on is not accounted for in this simulation. Calculation time: hours.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Protein data import; implemented by `options.pdb_mol=1`.
- Lines 19-20: Magnet field; implemented by `sys.magnet=21.1356`.
- Lines 22-23: Tolerances; implemented by `sys.tols.inter_cutoff=100`.
- Lines 26-27: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 38-39: Algorithmic options; implemented by `sys.enable={'prop_cache','op_cache','greedy'}`.
- Lines 41-42: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Kill carbons and nitrogens (protein assumed unlabelled); implemented by `spin_system=kill_spin(spin_system,strcmp('13C',spin_system.comp.isotopes))`.
- Lines 48-49: Build the basis; implemented by `spin_system=basis(spin_system,bas)`.
- Lines 51-52: Sequence parameters; implemented by `parameters.tmix=200e-3`.
- Lines 61-62: Simulation; implemented by `fid=liquid(spin_system,@noesy,parameters,'nmr')`.
- Lines 64-65: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 68-69: F2 Fourier transform; implemented by `f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.
- Lines 72-73: States signal; implemented by `f1_states=f1_cos-1i*f1_sin`.
- Lines 75-76: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 78-79: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 13: computes `options.pdb_mol` using `options.pdb_mol=1`.
- Lines 14: computes `options.select` using `options.select='all'`.
- Lines 15: computes `options.noshift` using `options.noshift='delete'`.
- Lines 16: computes `options.deuterate` using `options.deuterate='non-Me'`.
- Lines 17: computes `[sys,inter]` using `[sys,inter]=protein('2N9K.pdb','2N9K.bmrb',options)`.
- Lines 20: computes `sys.magnet` using `sys.magnet=21.1356`.
- Lines 23: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=100`.
- Lines 24: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=5.0`.
- Lines 27: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 28: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 29: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 30: computes `inter.tau_c` using `inter.tau_c={5e-9}`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 35: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 36: computes `bas.level` using `bas.level=2; bas.space_level=2`.
- Lines 39: computes `sys.enable` using `sys.enable={'prop_cache','op_cache','greedy'}`.
- Lines 42: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- 1H-1H NOESY spectrum of GB1 with everything deuterated except
- methyl groups. Deuteria are kept in the spin system because
- they are a part of the coupling network; methyl group rotati-
- on is not accounted for in this simulation.
- Calculation time: hours.
- Protein data import
- Magnet field
- Tolerances
- Relaxation theory
- Basis set
- Algorithmic options
- Create the spin system structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `protein()`, `create()`, `kill_spin()`, `strcmp()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
