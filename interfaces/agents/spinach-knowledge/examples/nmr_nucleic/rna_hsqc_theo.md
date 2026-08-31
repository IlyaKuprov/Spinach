# examples/nmr_nucleic/rna_hsqc_theo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_nucleic/rna_hsqc_theo.m`
- Signature: `rna_hsqc_theo()`
- Total lines: 79

## Purpose

1H-13C HSQC spectrum of the example RNA molecule provided by the Wagner group. Calculation time: minutes Shunsuke Imai Scott Robson Gerhard Wagner Zenawi Welderufael Ilya Kuprov

## Physical / mathematical content

- Nucleic-acid NMR examples. These files specialise biomolecular NMR workflows to RNA or DNA systems, with labelled nuclei, residue-level assignments, and multidimensional heteronuclear transfer logic.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Import RNA data; implemented by `options.noshift='delete'; options.deut_list={}`.
- Lines 18-19: Magnet field; implemented by `sys.magnet=17.62`.
- Lines 21-22: Tolerances; implemented by `sys.tols.inter_cutoff=5.0`.
- Lines 27-28: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 33-34: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 39-40: Sequence parameters; implemented by `parameters.J=90`.
- Lines 50-51: Create the spin system structure; implemented by `spin_system=create(sys,inter)`.
- Lines 53-54: Build the basis; implemented by `spin_system=basis(spin_system,bas)`.
- Lines 56-57: Simulation; implemented by `fid=liquid(spin_system,@hsqc,parameters,'nmr')`.
- Lines 59-60: Apodisation; implemented by `fid.pos=apodisation(spin_system,fid.pos,{{'cos'},{'cos'}})`.
- Lines 63-64: F2 Fourier transform; implemented by `f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1)`.
- Lines 67-68: Form States signal; implemented by `fid=f1_pos+conj(f1_neg)`.
- Lines 70-71: F1 Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2)`.
- Lines 73-74: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 15: computes `options.noshift` using `options.noshift='delete'; options.deut_list={}`.
- Lines 16: computes `[sys,inter]` using `[sys,inter]=nuclacid('example.pdb','example.txt',options)`.
- Lines 19: computes `sys.magnet` using `sys.magnet=17.62`.
- Lines 22: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=5.0`.
- Lines 23: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=5.0`.
- Lines 24: computes `sys.disable` using `sys.disable={'krylov','colorbar'}`.
- Lines 25: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 30: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 31: computes `bas.level` using `bas.level=4; bas.space_level=1`.
- Lines 34: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 35: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 36: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 37: computes `inter.damp_rate` using `inter.damp_rate=5.0`.
- Lines 40: computes `parameters.J` using `parameters.J=90`.
- Lines 41: computes `parameters.sweep` using `parameters.sweep=[2500 2000]`.
- Lines 42: computes `parameters.offset` using `parameters.offset=[22000 6000]`.

## Implementation structure

- 1H-13C HSQC spectrum of the example RNA molecule provided by
- the Wagner group.
- Calculation time: minutes
- Shunsuke Imai
- Scott Robson
- Gerhard Wagner
- Zenawi Welderufael
- Ilya Kuprov
- Import RNA data
- Magnet field
- Tolerances
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `nuclacid()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
