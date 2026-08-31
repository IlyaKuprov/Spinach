# examples/visualisation/hfc_pyrene.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/visualisation/hfc_pyrene.m`
- Signature: `hfc_pyrene()`
- Total lines: 26

## Purpose

Example of carbon hyperfine tensor visualisation for pyrene cation radical. Gaussian log is parsed.

## Physical / mathematical content

- Visualisation examples. These scripts expose tensor geometry: principal axes, shielding/hyperfine/EFG ellipsoids, molecular frames, and the relationship between tensor eigenstructure and observable anisotropy.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Read the Gaussian log; implemented by `props=gparse('pyrene_cation.log')`.
- Lines 11-12: Do the visualization; implemented by `kfigure()`.

### Key state/data transformations

- Lines 9: computes `props` using `props=gparse('pyrene_cation.log')`.
- Lines 14: computes `options.style` using `options.style='ellipsoids'`.

## Implementation structure

- Example of carbon hyperfine tensor visualisation for
- pyrene cation radical. Gaussian log is parsed.
- Read the Gaussian log
- Do the visualization

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gparse()`, `kfigure()`, `subplot()`, `hfc_display()`, `set()`, `ktitle()`, `scale_figure()`.
