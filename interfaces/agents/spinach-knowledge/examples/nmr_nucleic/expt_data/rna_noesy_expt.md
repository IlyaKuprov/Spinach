# examples/nmr_nucleic/expt_data/rna_noesy_expt.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_nucleic/expt_data/rna_noesy_expt.m`
- Signature: `rna_noesy_expt()`
- Total lines: 31

## Purpose

Experimental NOESY spectrum of the Harvard RNA. Shunsuke Imai Scott Robson Gerhard Wagner Zenawi Welderufael Ilya Kuprov

## Physical / mathematical content

- Nucleic-acid NMR examples. These files specialise biomolecular NMR workflows to RNA or DNA systems, with labelled nuclei, residue-level assignments, and multidimensional heteronuclear transfer logic.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnet field; implemented by `spin_system.inter.magnet=17.62`.
- Lines 14-15: Sequence parameters; implemented by `parameters.offset=3473`.
- Lines 21-22: Load data; implemented by `load('rna_noesy_expt.mat','spec_expt')`.
- Lines 24-25: Plotting; implemented by `spin_system.sys.disable={}; spin_system.sys.output=1`.

### Key state/data transformations

- Lines 12: computes `spin_system.inter.magnet` using `spin_system.inter.magnet=17.62`.
- Lines 15: computes `parameters.offset` using `parameters.offset=3473`.
- Lines 16: computes `parameters.sweep` using `parameters.sweep=[7500.00 7496.252]`.
- Lines 17: computes `parameters.zerofill` using `parameters.zerofill=[1024 4096]`.
- Lines 18: computes `parameters.spins` using `parameters.spins={'1H','1H'}`.
- Lines 19: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 25: computes `spin_system.sys.disable` using `spin_system.sys.disable={}; spin_system.sys.output=1`.

## Implementation structure

- Experimental NOESY spectrum of the Harvard RNA.
- Shunsuke Imai
- Scott Robson
- Gerhard Wagner
- Zenawi Welderufael
- Ilya Kuprov
- Magnet field
- Sequence parameters
- Load data
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
