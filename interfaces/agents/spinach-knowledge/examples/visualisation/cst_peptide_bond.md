# examples/visualisation/cst_peptide_bond.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/visualisation/cst_peptide_bond.m`
- Signature: `cst_peptide_bond()`
- Total lines: 33

## Purpose

Example of shielding tensor visualisation for a peptide bond. Gaussian log is parsed. Note: antisymmetric components of the shielding ten- sors are ignored.

## Physical / mathematical content

- Visualisation examples. These scripts expose tensor geometry: principal axes, shielding/hyperfine/EFG ellipsoids, molecular frames, and the relationship between tensor eigenstructure and observable anisotropy.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Read the Gaussian log; implemented by `props=gparse('../standard_systems/amino_acids/ala.log')`.
- Lines 14-15: Do the visualisation; implemented by `kfigure(); scale_figure([2.0 1.0])`.

### Key state/data transformations

- Lines 12: computes `props` using `props=gparse('../standard_systems/amino_acids/ala.log')`.
- Lines 17: computes `subplot(1,3,1); options.style` using `subplot(1,3,1); options.style='harmonics'`.
- Lines 22: computes `subplot(1,3,2); options.style` using `subplot(1,3,2); options.style='harmonics'`.
- Lines 27: computes `subplot(1,3,3); options.style` using `subplot(1,3,3); options.style='harmonics'`.

## Implementation structure

- Example of shielding tensor visualisation for a peptide
- bond. Gaussian log is parsed.
- Note: antisymmetric components of the shielding ten-
- sors are ignored.
- Read the Gaussian log
- Do the visualisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gparse()`, `kfigure()`, `scale_figure()`, `subplot()`, `cst_display()`, `set()`, `ktitle()`.
