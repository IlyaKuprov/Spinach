# examples/nmr_proteins/expt_data/noesy_ubiquitin_expt.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_proteins/expt_data/noesy_ubiquitin_expt.m`
- Signature: `noesy_ubiquitin_expt()`
- Total lines: 29

## Purpose

Experimental HNCO spectrum of human ubiquitin. Donghan Lee (Max Planck Institute) Ilya Kuprov (University of Southampton)

## Physical / mathematical content

- Protein NMR examples. These files specialise liquid-state pulse sequences to labelled biomolecules, exploiting one-bond and two-bond heteronuclear couplings, coherence pathway filtering, selective decoupling, and high-dimensional indirect detection.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Magnet field; implemented by `spin_system.inter.magnet=21.1356`.
- Lines 11-12: Sequence parameters; implemented by `parameters.offset=4250`.
- Lines 18-19: Load data; implemented by `load('noesy_ubiquitin_expt.mat','spectrum')`.
- Lines 21-22: Plotting; implemented by `spin_system.sys.disable={'colorbar'}`.

### Key state/data transformations

- Lines 9: computes `spin_system.inter.magnet` using `spin_system.inter.magnet=21.1356`.
- Lines 12: computes `parameters.offset` using `parameters.offset=4250`.
- Lines 13: computes `parameters.sweep` using `parameters.sweep=10815`.
- Lines 14: computes `parameters.zerofill` using `parameters.zerofill=[1024 1024]`.
- Lines 15: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 16: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 22: computes `spin_system.sys.disable` using `spin_system.sys.disable={'colorbar'}`.
- Lines 23: computes `spin_system.sys.output` using `spin_system.sys.output=1`.

## Implementation structure

- Experimental HNCO spectrum of human ubiquitin.
- Donghan Lee (Max Planck Institute)
- Ilya Kuprov (University of Southampton)
- Magnet field
- Sequence parameters
- Load data
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
