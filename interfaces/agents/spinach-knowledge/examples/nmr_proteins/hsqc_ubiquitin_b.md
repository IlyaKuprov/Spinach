# examples/nmr_proteins/hsqc_ubiquitin_b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/hsqc_ubiquitin_b.m`
- Signature: `hsqc_ubiquitin_b()`
- Total lines: 73

## Purpose

1H-15N HSQC of human ubiquitin, without 1H decoupling in F1 and 15N decoupling in F2. Nitrogen-proton multiplicity is retained in both dimensions. Calculation time: hours, faster with a Tesla A100 GPU. Zenawi Welderufael Luke Edwards Ilya Kuprov

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Protein data import; implemented by `options.pdb_mol=1`.
- Lines 19-20: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 22-23: Tolerances; implemented by `sys.tols.inter_cutoff=5.0`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 32-33: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Sequence parameters; implemented by `parameters.J=90`.
- Lines 50-51: Simulation; implemented by `fid=liquid(spin_system,@hsqc,parameters,'nmr')`.
- Lines 53-54: Apodisation; implemented by `fid.pos=apodisation(spin_system,fid.pos,{{'sqcos'},{'sqcos'}})`.
- Lines 57-58: F2 Fourier transform; implemented by `f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1)`.
- Lines 61-62: Form States signal; implemented by `fid=f1_pos+conj(f1_neg)`.
- Lines 64-65: F1 Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2)`.
- Lines 67-68: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 14: computes `options.pdb_mol` using `options.pdb_mol=1`.
- Lines 15: computes `options.noshift` using `options.noshift='delete'`.
- Lines 16: computes `options.select` using `options.select='backbone-hsqc'`.
- Lines 17: computes `[sys,inter]` using `[sys,inter]=protein('1D3Z.pdb','1D3Z.bmrb',options)`.
- Lines 20: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 23: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=5.0`.
- Lines 24: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 29: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 30: computes `bas.level` using `bas.level=4; bas.space_level=1`.
- Lines 33: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `parameters.J` using `parameters.J=90`.
- Lines 41: computes `parameters.sweep` using `parameters.sweep=[2400 4800]`.
- Lines 42: computes `parameters.offset` using `parameters.offset=[-7000 4500]`.
- Lines 43: computes `parameters.npoints` using `parameters.npoints=[256 256]`.
- Lines 44: computes `parameters.zerofill` using `parameters.zerofill=[1024 1024]`.

## Implementation structure

- 1H-15N HSQC of human ubiquitin, without 1H decoupling in F1 and
- 15N decoupling in F2. Nitrogen-proton multiplicity is retained
- in both dimensions.
- Calculation time: hours, faster with a Tesla A100 GPU.
- Zenawi Welderufael
- Luke Edwards
- Ilya Kuprov
- Protein data import
- Magnet field
- Tolerances
- Basis set
- Algorithmic options

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `protein()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
