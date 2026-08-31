# examples/nmr_proteins/hsqc_ubiquitin_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/hsqc_ubiquitin_a.m`
- Signature: `hsqc_ubiquitin_a()`
- Total lines: 72

## Purpose

1H-15N HSQC of human ubiquitin, decoupling applied in both dimensions. Calculation time: hours, faster with a Tesla A100 GPU. Zenawi Welderufael Luke Edwards Ilya Kuprov

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Protein data import; implemented by `options.pdb_mol=1`.
- Lines 18-19: Magnet field; implemented by `sys.magnet=11.7395`.
- Lines 21-22: Tolerances; implemented by `sys.tols.inter_cutoff=5.0`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Sequence parameters; implemented by `parameters.J=90`.
- Lines 49-50: Simulation; implemented by `fid=liquid(spin_system,@hsqc,parameters,'nmr')`.
- Lines 52-53: Apodisation; implemented by `fid.pos=apodisation(spin_system,fid.pos,{{'sqcos'},{'sqcos'}})`.
- Lines 56-57: F2 Fourier transform; implemented by `f1_pos=fftshift(fft(fid.pos,parameters.zerofill(2),1),1)`.
- Lines 60-61: Form States signal; implemented by `fid=f1_pos+conj(f1_neg)`.
- Lines 63-64: F1 Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill(1),2),2)`.
- Lines 66-67: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 13: computes `options.pdb_mol` using `options.pdb_mol=1`.
- Lines 14: computes `options.noshift` using `options.noshift='delete'`.
- Lines 15: computes `options.select` using `options.select='backbone-hsqc'`.
- Lines 16: computes `[sys,inter]` using `[sys,inter]=protein('1D3Z.pdb','1D3Z.bmrb',options)`.
- Lines 19: computes `sys.magnet` using `sys.magnet=11.7395`.
- Lines 22: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=5.0`.
- Lines 23: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 28: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 29: computes `bas.level` using `bas.level=4; bas.space_level=1`.
- Lines 32: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `parameters.J` using `parameters.J=90`.
- Lines 40: computes `parameters.sweep` using `parameters.sweep=[2000 4000]`.
- Lines 41: computes `parameters.offset` using `parameters.offset=[-5870 3753]`.
- Lines 42: computes `parameters.npoints` using `parameters.npoints=[128 256]`.
- Lines 43: computes `parameters.zerofill` using `parameters.zerofill=[1024 1024]`.

## Implementation structure

- 1H-15N HSQC of human ubiquitin, decoupling applied
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
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `protein()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
