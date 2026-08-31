# examples/nmr_nucleic/rna_noesy_theo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_nucleic/rna_noesy_theo.m`
- Signature: `rna_noesy_theo()`
- Total lines: 84

## Purpose

1H-1H NOESY spectrum of the example RNA molecule provided by the Gerhard Wagner group at Harvard University. Calculation time: hours Shunsuke Imai Scott Robson Gerhard Wagner Zenawi Welderufael Ilya Kuprov

## Physical / mathematical content

- Nucleic-acid NMR examples. These files specialise biomolecular NMR workflows to RNA or DNA systems, with labelled nuclei, residue-level assignments, and multidimensional heteronuclear transfer logic.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Import RNA data; implemented by `options.noshift='delete'`.
- Lines 20-21: Magnet field; implemented by `sys.magnet=17.62`.
- Lines 23-24: Tolerances; implemented by `sys.tols.inter_cutoff=1.0`.
- Lines 29-30: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 35-36: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 41-42: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Kill carbons and nitrogens (RNA assumed unlabelled); implemented by `spin_system=kill_spin(spin_system,strcmp('13C',spin_system.comp.isotopes))`.
- Lines 48-49: Build the basis; implemented by `spin_system=basis(spin_system,bas)`.
- Lines 51-52: Sequence parameters; implemented by `parameters.tmix=0.200`.
- Lines 61-62: Simulation; implemented by `fid=liquid(spin_system,@noesy,parameters,'nmr')`.
- Lines 64-65: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 68-69: F2 Fourier transform; implemented by `f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.
- Lines 72-73: States signal; implemented by `f1_states=f1_cos-1i*f1_sin`.
- Lines 75-76: F1 Fourier transform; implemented by `spec=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 78-79: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 15: computes `options.noshift` using `options.noshift='delete'`.
- Lines 16-17: computes `options.deut_list` using `options.deut_list={'GUA:H1','GUA:H21','GUA:H22','CYT:H41', 'CYT:H42','URI:H3','ADE:H61','ADE:H62'}`.
- Lines 18: computes `[sys,inter]` using `[sys,inter]=nuclacid('example.pdb','example.txt',options)`.
- Lines 21: computes `sys.magnet` using `sys.magnet=17.62`.
- Lines 24: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=1.0`.
- Lines 25: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=5.0`.
- Lines 26: computes `sys.disable` using `sys.disable={'krylov','colorbar'}`.
- Lines 27: computes `sys.enable` using `sys.enable={'prop_cache','greedy'}`.
- Lines 30: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 31: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 32: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 33: computes `inter.tau_c` using `inter.tau_c={3e-9}`.
- Lines 36: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 37: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 38: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 39: computes `bas.level` using `bas.level=5; bas.space_level=3`.
- Lines 42: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 52: computes `parameters.tmix` using `parameters.tmix=0.200`.

## Implementation structure

- 1H-1H NOESY spectrum of the example RNA molecule provided by
- the Gerhard Wagner group at Harvard University.
- Calculation time: hours
- Shunsuke Imai
- Scott Robson
- Gerhard Wagner
- Zenawi Welderufael
- Ilya Kuprov
- Import RNA data
- Magnet field
- Tolerances
- Relaxation theory

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `nuclacid()`, `create()`, `kill_spin()`, `strcmp()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
