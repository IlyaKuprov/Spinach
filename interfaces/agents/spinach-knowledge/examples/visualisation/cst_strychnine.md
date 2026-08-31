# examples/visualisation/cst_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/visualisation/cst_strychnine.m`
- Signature: `cst_strychnine()`
- Total lines: 28

## Purpose

Example of carbon shielding tensor visualisation for strychnine molecule. Gaussian log is parsed. Note: antisymmetric components of the shielding ten- sors are ignored.

## Physical / mathematical content

- Visualisation examples. These scripts expose tensor geometry: principal axes, shielding/hyperfine/EFG ellipsoids, molecular frames, and the relationship between tensor eigenstructure and observable anisotropy.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Read the Gaussian log; implemented by `props=gparse('strychnine.log')`.
- Lines 14-15: Do the visualisation; implemented by `kfigure(); subplot(1,2,1)`.

### Key state/data transformations

- Lines 12: computes `props` using `props=gparse('strychnine.log')`.
- Lines 16: computes `options.style` using `options.style='ellipsoids'`.

## Implementation structure

- Example of carbon shielding tensor visualisation for
- strychnine molecule. Gaussian log is parsed.
- Note: antisymmetric components of the shielding ten-
- sors are ignored.
- Read the Gaussian log
- Do the visualisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gparse()`, `kfigure()`, `subplot()`, `cst_display()`, `set()`, `ktitle()`, `scale_figure()`.
