# examples/nmr_proteins/hncoca_ubiquitin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/hncoca_ubiquitin.m`
- Signature: `hncoca_ubiquitin()`
- Total lines: 81

## Purpose

Theoretical HN(CO)CA of human ubiquitin. It is assumed that only the backbone is 13C,15N-labelled. Calculation time: minutes, faster with a Tesla A100 GPU.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Protein data import; implemented by `options.pdb_mol=1`.
- Lines 17-18: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 20-21: Tolerances; implemented by `sys.tols.inter_cutoff=2.0`.
- Lines 24-25: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Sequence parameters; implemented by `parameters.tau=[2.25e-3, 2.75e-3, 8.00e-3, 7.00e-3]`.
- Lines 47-48: Simulation; implemented by `fid=liquid(spin_system,@hncoca,parameters,'nmr')`.
- Lines 50-51: Apodisation; implemented by `fid.pos_pos=apodisation(spin_system,fid.pos_pos,{{'sqcos'},{'sqcos'},{'sqcos'}})`.
- Lines 56-57: F3 Fourier transform; implemented by `f3_pos_pos=fftshift(fft(fid.pos_pos,parameters.zerofill(3),3),3)`.
- Lines 62-63: Absorption part of F3 signal; implemented by `f3_pos=f3_pos_pos+conj(f3_neg_neg)`.
- Lines 66-67: F2 Fourier transform; implemented by `f3f2_pos=fftshift(fft(f3_pos,parameters.zerofill(2),2),2)`.
- Lines 70-71: Absorption part of F2 signal; implemented by `f3f2=f3f2_pos+conj(f3f2_neg)`.
- Lines 73-74: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f3f2,parameters.zerofill(1),1),1)`.
- Lines 76-78: Plotting; implemented by `kfigure(); plot_3d(spin_system,-real(spectrum),parameters, 10,[0.2 0.9 0.2 0.9],2,'positive')`.

### Key state/data transformations

- Lines 12: computes `options.pdb_mol` using `options.pdb_mol=1`.
- Lines 13: computes `options.noshift` using `options.noshift='delete'`.
- Lines 14: computes `options.select` using `options.select='backbone-minimal'`.
- Lines 15: computes `[sys,inter]` using `[sys,inter]=protein('1D3Z.pdb','1D3Z.bmrb',options)`.
- Lines 18: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 21: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=2.0`.
- Lines 22: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 25: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 27: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 28: computes `bas.level` using `bas.level=4; bas.space_level=1`.
- Lines 31: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 32: computes `sys.disable` using `sys.disable={'krylov'}`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `parameters.tau` using `parameters.tau=[2.25e-3, 2.75e-3, 8.00e-3, 7.00e-3]`.
- Lines 40: computes `parameters.spins` using `parameters.spins={'15N','13C','1H'}`.
- Lines 41: computes `parameters.offset` using `parameters.offset=[-7100 8450 4850]`.
- Lines 42: computes `parameters.sweep` using `parameters.sweep=[2500 4500 3000]`.

## Implementation structure

- Theoretical HN(CO)CA of human ubiquitin. It is assumed that
- only the backbone is 13C,15N-labelled.
- Calculation time: minutes, faster with a Tesla A100 GPU.
- Protein data import
- Magnet field
- Tolerances
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `protein()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `plot_3d()`.
