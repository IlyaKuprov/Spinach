# examples/nmr_proteins/ct_cosy_gb1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/ct_cosy_gb1.m`
- Signature: `ct_cosy_gb1()`
- Total lines: 67

## Purpose

Constant-time COSY experiment simulation for the GB1 protein. Simulation time: minutes, faster with a Tesla A100 GPU.

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Protein data import; implemented by `options.pdb_mol=1`.
- Lines 17-18: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 20-21: Tolerances; implemented by `sys.tols.inter_cutoff=1.0`.
- Lines 24-25: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 33-34: Sequence parameters; implemented by `parameters.offset=2400`.
- Lines 41-42: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Kill carbons and nitrogens (protein assumed unlabelled); implemented by `spin_system=kill_spin(spin_system,strcmp('13C',spin_system.comp.isotopes))`.
- Lines 48-49: Build the basis; implemented by `spin_system=basis(spin_system,bas)`.
- Lines 51-52: Simulation; implemented by `fid=liquid(spin_system,@ct_cosy,parameters,'nmr')`.
- Lines 54-55: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sqcos'},{'sqcos'}})`.
- Lines 57-58: Fourier transform; implemented by `spectrum=fftn(fid,parameters.zerofill)`.
- Lines 61-62: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 12: computes `options.pdb_mol` using `options.pdb_mol=1`.
- Lines 13: computes `options.noshift` using `options.noshift='delete'`.
- Lines 14: computes `options.select` using `options.select='backbone'`.
- Lines 15: computes `[sys,inter]` using `[sys,inter]=protein('2N9K.pdb','2N9K.bmrb',options)`.
- Lines 18: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 21: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=1.0`.
- Lines 22: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=3.0`.
- Lines 25: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 27: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 28: computes `bas.level` using `bas.level=4; bas.space_level=1`.
- Lines 31: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 34: computes `parameters.offset` using `parameters.offset=2400`.
- Lines 35: computes `parameters.sweep` using `parameters.sweep=[9000 9000]`.
- Lines 36: computes `parameters.npoints` using `parameters.npoints=[256 256]`.
- Lines 37: computes `parameters.zerofill` using `parameters.zerofill=[512 512]`.
- Lines 38: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 39: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.

## Implementation structure

- Constant-time COSY experiment simulation for the
- GB1 protein.
- Simulation time: minutes, faster with a Tesla A100 GPU.
- Protein data import
- Magnet field
- Tolerances
- Basis set
- Algorithmic options
- Sequence parameters
- Spinach housekeeping
- Kill carbons and nitrogens (protein assumed unlabelled)
- Build the basis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `protein()`, `create()`, `kill_spin()`, `strcmp()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
